# Sequence Alignment Algorithms Comparison

This project compares two sequence alignment algorithms: a naive approach and a suffix array-based approach.


## Description

- `naive.py`: Implements a simple, brute-force alignment algorithm.
- `star.py`: Implements a suffix array-based alignment algorithm.
- `benchmark.py`: Runs performance tests on both algorithms.
- `generate_fasta.py`: Generates test data (reference and reads).


## Running the Benchmark

To run the benchmark:

```
python scripts/benchmark.py
```

This will execute both algorithms multiple times and provide average execution times and standard deviations.

## Generating Test Data

To generate new test data:

```
python scripts/generate_fasta.py
```

This will create new reference and reads files in the `data/` directory.
