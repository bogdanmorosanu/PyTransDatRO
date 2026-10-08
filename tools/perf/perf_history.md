# PyTransDatRO Performance & Optimization History Log
This document tracks performance metrics over time across versions and optimization experiments.
*Last Updated:* `2026-10-08 21:13:33`

## 1. Run Registry

| Run ID | Version / Tag | Timestamp | Description | Init (ms) |
|---|---|---|---|---|
| `run_1791485258` | **baseline** | 2026-10-08T19:47:38 | Initial unoptimized baseline | 6.54 |
| `run_1791487611` | **opt_phase1_2** | 2026-10-08T20:26:51 | Implemented bulk processing (_many methods) + minor inline math opts | 6.08 |
| `run_1791487666` | **opt_phase3** | 2026-10-08T20:27:46 | Increased LRU cache sizes to cover the entire grid space | 5.86 |
| `run_1791489344` | **opt_polymorphic** | 2026-10-08T20:55:44 | Unified polymorphic API (scalars and sequences) with removed _many duplicates | 8.60 |
| `run_1791490413` | **opt_approach_2** | 2026-10-08T21:13:33 | Refactored with Approach 2: Boundary normalization in TransRO, strictly single-implementation sequence methods internally | 5.68 |

---
## 2. Latest Benchmark Summary (`opt_approach_2`)

**Description:** Refactored with Approach 2: Boundary normalization in TransRO, strictly single-implementation sequence methods internally  
**Python:** 3.13.16 (Windows-11-10.0.26300-SP0)  
**Cold Initialization:** `5.68 ms`

### Detailed Scenario Metrics

| Scenario | Points | Direction | Mode | Cache | Best Time (s) | Throughput (pts/s) | Latency (µs/pt) |
|---|---|---|---|---|---|---|---|
| S1_micro | 10 | st70_to_etrs89 | 3D | Warm | `0.0002` | **47,237** | 21.2 |
| S1_micro | 10 | etrs89_to_st70 | 3D | Warm | `0.0002` | **50,891** | 19.6 |
| S2_small | 250 | st70_to_etrs89 | 3D | Cold | `0.0051` | **48,625** | 20.6 |
| S2_small | 250 | st70_to_etrs89 | 3D | Warm | `0.0050` | **49,818** | 20.1 |
| S2_small | 250 | etrs89_to_st70 | 3D | Cold | `0.0044` | **56,325** | 17.8 |
| S2_small | 250 | etrs89_to_st70 | 3D | Warm | `0.0044` | **56,414** | 17.7 |
| S3_medium | 2500 | st70_to_etrs89 | 3D | Cold | `0.0509` | **49,099** | 20.4 |
| S3_medium | 2500 | st70_to_etrs89 | 3D | Warm | `0.0501` | **49,908** | 20.0 |
| S3_medium | 2500 | etrs89_to_st70 | 3D | Cold | `0.0466` | **53,699** | 18.6 |
| S3_medium | 2500 | etrs89_to_st70 | 3D | Warm | `0.0457` | **54,714** | 18.3 |
| S4_large | 25000 | st70_to_etrs89 | 2D | Warm | `0.5800` | **43,104** | 23.2 |
| S4_large | 25000 | st70_to_etrs89 | 3D | Cold | `0.5933` | **42,137** | 23.7 |
| S4_large | 25000 | st70_to_etrs89 | 3D | Warm | `0.5554` | **45,010** | 22.2 |
| S4_large | 25000 | etrs89_to_st70 | 2D | Warm | `0.5069` | **49,322** | 20.3 |
| S4_large | 25000 | etrs89_to_st70 | 3D | Cold | `0.5764` | **43,373** | 23.1 |
| S4_large | 25000 | etrs89_to_st70 | 3D | Warm | `0.5040` | **49,606** | 20.2 |

---
## 3. Comparative Evolution vs Baseline (`baseline`)

| Scenario & Mode | Direction | Cache | Baseline Throughput | Latest Throughput | Speedup (Delta %) |
|---|---|---|---|---|---|
| S1_micro (3D) | st70_to_etrs89 | Warm | 39,920 pts/s | 47,237 pts/s | **+18.3%** |
| S1_micro (3D) | etrs89_to_st70 | Warm | 46,642 pts/s | 50,891 pts/s | **+9.1%** |
| S2_small (3D) | st70_to_etrs89 | Cold | 39,976 pts/s | 48,625 pts/s | **+21.6%** |
| S2_small (3D) | st70_to_etrs89 | Warm | 42,257 pts/s | 49,818 pts/s | **+17.9%** |
| S2_small (3D) | etrs89_to_st70 | Cold | 48,535 pts/s | 56,325 pts/s | **+16.1%** |
| S2_small (3D) | etrs89_to_st70 | Warm | 49,350 pts/s | 56,414 pts/s | **+14.3%** |
| S3_medium (3D) | st70_to_etrs89 | Cold | 40,218 pts/s | 49,099 pts/s | **+22.1%** |
| S3_medium (3D) | st70_to_etrs89 | Warm | 41,291 pts/s | 49,908 pts/s | **+20.9%** |
| S3_medium (3D) | etrs89_to_st70 | Cold | 46,780 pts/s | 53,699 pts/s | **+14.8%** |
| S3_medium (3D) | etrs89_to_st70 | Warm | 47,549 pts/s | 54,714 pts/s | **+15.1%** |
| S4_large (2D) | st70_to_etrs89 | Warm | 27,308 pts/s | 43,104 pts/s | **+57.8%** |
| S4_large (3D) | st70_to_etrs89 | Cold | 26,790 pts/s | 42,137 pts/s | **+57.3%** |
| S4_large (3D) | st70_to_etrs89 | Warm | 25,673 pts/s | 45,010 pts/s | **+75.3%** |
| S4_large (2D) | etrs89_to_st70 | Warm | 31,896 pts/s | 49,322 pts/s | **+54.6%** |
| S4_large (3D) | etrs89_to_st70 | Cold | 30,265 pts/s | 43,373 pts/s | **+43.3%** |
| S4_large (3D) | etrs89_to_st70 | Warm | 33,340 pts/s | 49,606 pts/s | **+48.8%** |
