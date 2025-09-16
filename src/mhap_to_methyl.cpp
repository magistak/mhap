#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <cinttypes>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <zlib.h>

namespace {

struct StringHash {
    using is_transparent = void;
    std::size_t operator()(std::string_view s) const noexcept {
        std::size_t h = 1469598103934665603ULL;
        for (unsigned char c : s) {
            h ^= c;
            h *= 1099511628211ULL;
        }
        return h;
    }
    std::size_t operator()(const std::string &s) const noexcept { return operator()(std::string_view{s}); }
    std::size_t operator()(const char *s) const noexcept { return operator()(std::string_view{s}); }
};

struct StringEq {
    using is_transparent = void;
    bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
    bool operator()(const std::string &a, const std::string &b) const noexcept { return a == b; }
    bool operator()(const char *a, const std::string &b) const noexcept { return b == a; }
    bool operator()(const std::string &a, const char *b) const noexcept { return a == b; }
};

struct ChromosomeStats {
    std::string name;
    std::vector<uint32_t> positions;
};

struct GenomeData {
    std::vector<ChromosomeStats> chromosomes;
    std::unordered_map<std::string, std::size_t, StringHash, StringEq> lookup;
};

struct WarningCounters {
    std::atomic<uint64_t> count_le_zero{0};
    std::atomic<uint64_t> unknown_chrom{0};
    std::atomic<uint64_t> interval_no_cpg{0};
    std::atomic<uint64_t> hap_mismatch{0};
    std::atomic<uint64_t> invalid_hap_char{0};
};

class GzReader {
public:
    explicit GzReader(const std::string &path) {
        file_ = gzopen(path.c_str(), "rb");
        if (!file_) {
            throw std::runtime_error("Failed to open gzip file: " + path);
        }
    }

    ~GzReader() {
        if (file_) {
            gzclose(file_);
        }
    }

    bool getline(std::string &out) {
        out.clear();
        const int buffer_size = 1 << 14;
        char buffer[buffer_size];
        bool read_any = false;
        while (true) {
            char *res = gzgets(file_, buffer, buffer_size);
            if (!res) {
                if (!read_any) {
                    return false;
                }
                return true;
            }
            read_any = true;
            int len = static_cast<int>(std::strlen(buffer));
            bool ended = false;
            if (len > 0 && buffer[len - 1] == '\n') {
                ended = true;
                --len;
            }
            if (len > 0 && buffer[len - 1] == '\r') {
                --len;
            }
            out.append(buffer, len);
            if (ended) {
                return true;
            }
        }
    }

private:
    gzFile file_;
};

bool parse_uint32(std::string_view value, uint32_t &out) {
    if (value.empty()) {
        return false;
    }
    uint64_t acc = 0;
    for (char ch : value) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        acc = acc * 10 + static_cast<unsigned>(ch - '0');
        if (acc > std::numeric_limits<uint32_t>::max()) {
            return false;
        }
    }
    out = static_cast<uint32_t>(acc);
    return true;
}

bool parse_int64(std::string_view value, int64_t &out) {
    if (value.empty()) {
        return false;
    }
    bool neg = false;
    std::size_t idx = 0;
    if (value[0] == '-') {
        neg = true;
        idx = 1;
        if (idx == value.size()) {
            return false;
        }
    }
    int64_t acc = 0;
    for (; idx < value.size(); ++idx) {
        char ch = value[idx];
        if (ch < '0' || ch > '9') {
            return false;
        }
        int digit = ch - '0';
        if (neg) {
            if (acc < (std::numeric_limits<int64_t>::lowest() + digit) / 10) {
                return false;
            }
            acc = acc * 10 - digit;
        } else {
            if (acc > (std::numeric_limits<int64_t>::max() - digit) / 10) {
                return false;
            }
            acc = acc * 10 + digit;
        }
    }
    out = acc;
    return true;
}

bool split_fields(const std::string &line, std::array<std::string_view, 6> &fields) {
    std::size_t start = 0;
    for (int i = 0; i < 5; ++i) {
        std::size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            return false;
        }
        fields[i] = std::string_view(line.data() + start, tab - start);
        start = tab + 1;
    }
    fields[5] = std::string_view(line.data() + start, line.size() - start);
    return true;
}

GenomeData load_cpg_index(const std::string &path) {
    GenomeData data;
    GzReader reader(path);
    std::string line;
    ChromosomeStats *current = nullptr;
    while (reader.getline(line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::size_t tab = line.find('\t');
        if (tab == std::string::npos) {
            throw std::runtime_error("Malformed CpG index line: " + line);
        }
        std::size_t tab2 = line.find('\t', tab + 1);
        std::string_view chrom(line.data(), tab);
        std::size_t pos_len = (tab2 == std::string::npos ? line.size() : tab2) - (tab + 1);
        std::string_view pos_str(line.data() + tab + 1, pos_len);
        uint32_t pos = 0;
        if (!parse_uint32(pos_str, pos)) {
            throw std::runtime_error("Invalid CpG position: " + std::string(pos_str));
        }
        if (!current || current->name != chrom) {
            data.chromosomes.emplace_back();
            ChromosomeStats &chr = data.chromosomes.back();
            chr.name.assign(chrom.begin(), chrom.end());
            data.lookup.emplace(chr.name, data.chromosomes.size() - 1);
            current = &chr;
        }
        current->positions.push_back(pos);
    }
    return data;
}

struct CLIOptions {
    std::string cpg_path;
    std::string output_path;
    std::string output_dir;
    std::string output_suffix = ".tsv.gz";
    std::vector<std::string> mhap_paths;
    std::size_t threads = 1;
    uint64_t min_cov = 1;
};

bool ends_with(const std::string &value, const char *suffix) {
    std::size_t len = std::strlen(suffix);
    if (value.size() < len) {
        return false;
    }
    return std::equal(suffix, suffix + len, value.end() - len);
}

bool strip_suffix(std::string &value, std::string_view suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    if (value.compare(value.size() - suffix.size(), suffix.size(), suffix.data(), suffix.size()) == 0) {
        value.resize(value.size() - suffix.size());
        return true;
    }
    return false;
}

std::string basename(std::string_view path) {
    std::size_t pos = path.find_last_of("/\\");
    if (pos == std::string_view::npos) {
        return std::string(path);
    }
    return std::string(path.substr(pos + 1));
}

void append_mhap_list(const std::string &list_path, std::vector<std::string> &out) {
    std::ifstream in(list_path);
    if (!in) {
        throw std::runtime_error("Failed to open mhap list: " + list_path);
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) {
            continue;
        }
        std::size_t end = line.find_last_not_of(" \t\r\n");
        out.emplace_back(line.substr(start, end - start + 1));
    }
}

void print_usage(const char *prog) {
    std::cerr << "Usage: " << prog << " --cpg <hg19_CpG.gz> [--output <file>|--output-dir <dir>] [--output-suffix .ext] [--threads N] [--min-cov C] [--mhap-list list.txt] <mhap1.gz> [mhap2.gz ...]\n";
}

std::optional<CLIOptions> parse_cli(int argc, char **argv) {
    CLIOptions opts;
    if (argc < 5) {
        print_usage(argv[0]);
        return std::nullopt;
    }
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--cpg" && i + 1 < argc) {
            opts.cpg_path = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            opts.output_path = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            opts.output_dir = argv[++i];
        } else if (arg == "--output-suffix" && i + 1 < argc) {
            opts.output_suffix = argv[++i];
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            int64_t val = 0;
            if (!parse_int64(argv[i + 1], val) || val <= 0) {
                std::cerr << "Invalid thread count: " << argv[i + 1] << "\n";
                return std::nullopt;
            }
            opts.threads = static_cast<std::size_t>(val);
            ++i;
        } else if (arg == "--min-cov" && i + 1 < argc) {
            int64_t val = 0;
            if (!parse_int64(argv[i + 1], val) || val < 0) {
                std::cerr << "Invalid min coverage: " << argv[i + 1] << "\n";
                return std::nullopt;
            }
            opts.min_cov = static_cast<uint64_t>(val);
            ++i;
        } else if (arg == "--mhap-list" && i + 1 < argc) {
            append_mhap_list(argv[++i], opts.mhap_paths);
        } else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return std::nullopt;
        } else {
            opts.mhap_paths.emplace_back(std::move(arg));
        }
    }
    if (!opts.output_path.empty() && !opts.output_dir.empty()) {
        std::cerr << "Specify either --output or --output-dir, not both.\n";
        return std::nullopt;
    }
    if (opts.output_path.empty() && opts.output_dir.empty()) {
        std::cerr << "Missing --output or --output-dir destination.\n";
        return std::nullopt;
    }
    if (opts.cpg_path.empty() || opts.mhap_paths.empty()) {
        print_usage(argv[0]);
        return std::nullopt;
    }
    if (opts.output_suffix.empty()) {
        std::cerr << "Output suffix must be non-empty.\n";
        return std::nullopt;
    }
    if (opts.mhap_paths.size() > 1 && opts.output_dir.empty()) {
        std::cerr << "Multiple input files require --output-dir.\n";
        return std::nullopt;
    }
    if (opts.threads == 0) {
        opts.threads = 1;
    }
    return opts;
}

std::string derive_output_path(const CLIOptions &opts, const std::string &input_path) {
    if (!opts.output_dir.empty()) {
        std::string name = basename(input_path);
        strip_suffix(name, ".bgz");
        strip_suffix(name, ".gz");
        strip_suffix(name, ".mhap");
        if (name.empty()) {
            name = "sample";
        }
        std::string path = opts.output_dir;
        if (!path.empty() && path.back() != '/' && path.back() != '\\') {
            path.push_back('/');
        }
        path += name;
        path += opts.output_suffix;
        return path;
    }
    return opts.output_path;
}

class OutputWriter {
public:
    explicit OutputWriter(const std::string &path) : path_(path) {
        if (ends_with(path, ".gz")) {
            gz_ = gzopen(path.c_str(), "wb");
            if (!gz_) {
                throw std::runtime_error("Failed to open gzip output: " + path);
            }
        } else {
            file_.open(path);
            if (!file_) {
                throw std::runtime_error("Failed to open output file: " + path);
            }
            file_ << std::fixed << std::setprecision(6);
        }
    }

    ~OutputWriter() {
        if (gz_) {
            gzclose(gz_);
        }
    }

    void write_line(std::string_view line) { write_bytes(line.data(), line.size()); }

    void write_bytes(const char *data, std::size_t len) {
        if (gz_) {
            if (gzwrite(gz_, data, static_cast<unsigned>(len)) != static_cast<int>(len)) {
                throw std::runtime_error("Failed to write gzip output: " + path_);
            }
        } else {
            file_.write(data, static_cast<std::streamsize>(len));
        }
    }

private:
    std::string path_;
    gzFile gz_{nullptr};
    std::ofstream file_;
};

void write_output(const GenomeData &genome, const std::vector<std::vector<uint64_t>> &coverage,
                  const std::vector<std::vector<uint64_t>> &methyl, const std::string &output_path,
                  uint64_t min_cov) {
    OutputWriter writer(output_path);
    writer.write_line("chrom\tpos\tcov\tmeth\tbeta\n");
    char buffer[128];
    for (std::size_t chr_idx = 0; chr_idx < genome.chromosomes.size(); ++chr_idx) {
        const ChromosomeStats &chr = genome.chromosomes[chr_idx];
        const auto &cov_vec = coverage[chr_idx];
        const auto &meth_vec = methyl[chr_idx];
        for (std::size_t i = 0; i < chr.positions.size(); ++i) {
            uint64_t cov = cov_vec[i];
            if (cov < min_cov) {
                continue;
            }
            uint64_t meth = meth_vec[i];
            double beta = cov == 0 ? 0.0 : static_cast<double>(meth) / static_cast<double>(cov);
            int n = std::snprintf(buffer, sizeof(buffer), "%s\t%u\t%" PRIu64 "\t%" PRIu64 "\t%.6f\n", chr.name.c_str(), chr.positions[i], cov, meth, beta);
            if (n < 0 || static_cast<std::size_t>(n) >= sizeof(buffer)) {
                throw std::runtime_error("Formatting error while writing output");
            }
            writer.write_bytes(buffer, static_cast<std::size_t>(n));
        }
    }
}

void process_file(const std::string &path, const GenomeData &genome, WarningCounters &warnings,
                  const CLIOptions &opts) {
    GzReader reader(path);
    std::string line;
    std::vector<std::vector<uint64_t>> coverage(genome.chromosomes.size());
    std::vector<std::vector<uint64_t>> methyl(genome.chromosomes.size());
    for (std::size_t i = 0; i < genome.chromosomes.size(); ++i) {
        std::size_t n = genome.chromosomes[i].positions.size();
        coverage[i].assign(n, 0);
        methyl[i].assign(n, 0);
    }
    while (reader.getline(line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::array<std::string_view, 6> fields;
        if (!split_fields(line, fields)) {
            continue;
        }
        std::string_view chrom = fields[0];
        int64_t start = 0;
        int64_t end = 0;
        int64_t count = 0;
        if (!parse_int64(fields[1], start) || !parse_int64(fields[2], end) || !parse_int64(fields[4], count)) {
            continue;
        }
        if (count <= 0) {
            warnings.count_le_zero.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        std::string chrom_key(chrom);
        auto it = genome.lookup.find(chrom_key);
        if (it == genome.lookup.end()) {
            warnings.unknown_chrom.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        std::size_t chr_index = it->second;
        const auto &positions = genome.chromosomes[chr_index].positions;
        if (positions.empty()) {
            continue;
        }
        if (end < positions.front() || start > positions.back()) {
            warnings.interval_no_cpg.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        if (start < 1 || end < start) {
            continue;
        }
        auto lower = std::lower_bound(positions.begin(), positions.end(), static_cast<uint32_t>(start));
        if (lower == positions.end()) {
            warnings.interval_no_cpg.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        auto upper = std::upper_bound(lower, positions.end(), static_cast<uint32_t>(end));
        std::size_t idx_start = static_cast<std::size_t>(lower - positions.begin());
        std::size_t idx_end = static_cast<std::size_t>(upper - positions.begin());
        if (idx_end <= idx_start) {
            warnings.interval_no_cpg.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        std::size_t expected = idx_end - idx_start;
        std::string_view haplotype = fields[3];
        if (haplotype.size() != expected) {
            warnings.hap_mismatch.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        auto &cov_vec = coverage[chr_index];
        auto &meth_vec = methyl[chr_index];
        uint64_t count_u = static_cast<uint64_t>(count);
        for (std::size_t offset = 0; offset < expected; ++offset) {
            char base = haplotype[offset];
            std::size_t idx = idx_start + offset;
            cov_vec[idx] += count_u;
            if (base == '1') {
                meth_vec[idx] += count_u;
            } else if (base == '0') {
                // unmethylated
            } else {
                warnings.invalid_hap_char.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
    std::string output_path = derive_output_path(opts, path);
    write_output(genome, coverage, methyl, output_path, opts.min_cov);
}

void print_warnings(const WarningCounters &warnings) {
    auto print_if = [](const char *label, uint64_t value) {
        if (value) {
            std::cerr << label << ": " << value << "\n";
        }
    };
    print_if("Records_skipped_count_le_zero", warnings.count_le_zero.load());
    print_if("Records_skipped_unknown_chrom", warnings.unknown_chrom.load());
    print_if("Intervals_without_CpG", warnings.interval_no_cpg.load());
    print_if("Haplotype_length_mismatches", warnings.hap_mismatch.load());
    print_if("Invalid_haplotype_characters", warnings.invalid_hap_char.load());
}

} // namespace

int main(int argc, char **argv) {
    auto options = parse_cli(argc, argv);
    if (!options) {
        return 1;
    }
    CLIOptions opts = std::move(*options);
    try {
        GenomeData genome = load_cpg_index(opts.cpg_path);
        if (genome.chromosomes.empty()) {
            std::cerr << "No CpG entries found." << std::endl;
            return 1;
        }
        WarningCounters warnings;
        std::vector<std::thread> workers;
        std::mutex queue_mutex;
        std::size_t next_index = 0;
        std::size_t thread_count = std::max<std::size_t>(1, opts.threads);
        workers.reserve(thread_count);
        auto worker = [&]() {
            while (true) {
                std::size_t idx = 0;
                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    if (next_index >= opts.mhap_paths.size()) {
                        break;
                    }
                    idx = next_index++;
                }
                try {
                    process_file(opts.mhap_paths[idx], genome, warnings, opts);
                } catch (const std::exception &ex) {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    std::cerr << "Error processing " << opts.mhap_paths[idx] << ": " << ex.what() << "\n";
                }
            }
        };
        for (std::size_t i = 0; i < thread_count; ++i) {
            workers.emplace_back(worker);
        }
        for (auto &t : workers) {
            t.join();
        }
        print_warnings(warnings);
    } catch (const std::exception &ex) {
        std::cerr << "Fatal error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
