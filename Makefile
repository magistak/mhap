CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -march=native -Wall -Wextra -Wpedantic
LDFLAGS ?=
LIBS = -lz -pthread

BIN = mhap_to_methyl
SRC = src/mhap_to_methyl.cpp

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -f $(BIN)

.PHONY: all clean
