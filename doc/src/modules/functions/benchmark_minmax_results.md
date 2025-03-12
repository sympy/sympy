# Min/Max Performance with evaluate=False

This document provides benchmark results showing the performance improvement when using `evaluate=False` with Min/Max functions that have many arguments.

## Background

Issue #16249 reported that Min/Max with many arguments was slow to construct, and setting `evaluate=False` didn't help because it didn't disable the `_find_localzeros` function.

The issue has been fixed in the current implementation, where `_find_localzeros` is now properly skipped when `evaluate=False` is set.

## Benchmark Results

The following benchmark was run with 50 symbols:

Running benchmark with 50 symbols, 3 runs each...
Testing Min with evaluate=True...
Run 1: 1.4782 seconds
Run 2: 1.1790 seconds
Run 3: 1.2110 seconds
Average time with Min(evaluate=True): 1.2894 seconds
Testing Min with evaluate=False...
Run 1: 0.0010 seconds
Run 2: 0.0000 seconds
Run 3: 0.0000 seconds
Average time with Min(evaluate=False): 0.0003 seconds
Testing Max with evaluate=True...
Run 1: 1.1770 seconds
Run 2: 1.1460 seconds
Run 3: 1.1250 seconds
Average time with Max(evaluate=True): 1.1493 seconds
Testing Max with evaluate=False...
Run 1: 0.0010 seconds
Run 2: 0.0010 seconds
Run 3: 0.0000 seconds
Average time with Max(evaluate=False): 0.0007 seconds
=== SUMMARY ===
Min with evaluate=True:  1.2894 seconds
Min with evaluate=False: 0.0003 seconds
Min speedup factor: 3870.45x
Max with evaluate=True:  1.1493 seconds
Max with evaluate=False: 0.0007 seconds
Max speedup factor: 1727.22x

## Conclusion

The benchmark results show that using `evaluate=False` with Min/Max functions now provides a significant performance improvement (over 1000x speedup) when dealing with many arguments.

This confirms that the issue #16249 has been fixed in the current implementation.