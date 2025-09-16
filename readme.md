**Inputs**

1. CpG index file (hg19_CpG.gz)
    
    - Tabix-indexed, bgzip-compressed.
    - Each data line: chrom<TAB>pos
        - chrom: chromosome name (e.g. chr1)
        - pos: 1-based genomic coordinate of a CpG cytosine (the C of a CpG on the + strand).
    - Must be sorted by chrom then pos.
    - Provides the complete ordered list of CpG positions per chromosome.
2. One or more mHap files (*.mhap.gz)
    
    - Bgzip-compressed text.
    - Each non-comment line (6 columns):
        1. chrom (string)
        2. start (1-based inclusive integer)
        3. end (1-based inclusive integer)
        4. haplotype (string of characters '0' and '1'; length = number of CpG sites represented)
        5. count (integer ≥ 1; number of molecules/reads supporting this haplotype)
        6. strand ('+' or '-')
    - Interpretation:
        - The haplotype string encodes ONLY the CpG sites that lie within [start, end].
        - Position i in haplotype corresponds to the i-th CpG (in ascending genomic order) whose coordinate is between start and end inclusive.
        - '1' = methylated, '0' = unmethylated.
        - Length of haplotype should equal the number of CpGs in the genomic interval; if it differs, there should be a warning

**Output Columns (row per CpG site with coverage ≥ min_cov threshold)**

- chrom: chromosome name
- pos: 1-based CpG genomic coordinate
- cov: total supporting count (sum of counts of all haplotypes covering this CpG)
- meth: total methylated count (subset of cov where haplotype char was '1')
- beta: meth / cov (floating point; 0 if cov == 0)

**Edge / Data Quality Notes**

- Records with count ≤ 0 are ignored.
- Non '0'/'1' characters (if any) should be treated: warning and treat as 0
- If haplotype length ≠ number of CpGs in interval: warning