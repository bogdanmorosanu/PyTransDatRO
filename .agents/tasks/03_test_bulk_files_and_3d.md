# Bulk Transformation and 3D Test Implementation

## Objectives
- [x] Refactor existing Stereo70 $\rightarrow$ ETRS89 bulk test to use parsed session-scoped fixtures.
- [x] Implement inverse bulk test (ETRS89 $\rightarrow$ Stereo70).
- [x] Implement 3D (Elevation) test for Stereo70 $\rightarrow$ ETRS89 (marked as xfail due to known grid discrepancies).
- [x] Implement 3D (Elevation) test for ETRS89 $\rightarrow$ Stereo70 (marked as xfail due to known grid discrepancies).

## Implementation Details
1. Created new fixtures (`st70_to_etrs89_input_data`, `st70_to_etrs89_expected_data`, `etrs89_to_st70_input_data`, `etrs89_to_st70_expected_data`) in `tests/test_trans_ro.py` using `@pytest.fixture(scope="session")` to load and parse CSV data directly into lists of tuples. This optimization prevents redundant file I/O operations and casting during each test.
2. Updated bulk testing functions to use the session data and verified they work correctly for planimetric (2D) transformations.
3. Added corresponding 3D tests and configured them with `@pytest.mark.xfail` based on the user's notification of expected grid elevation mismatches at this stage.

## Completion Status
The tests are running smoothly and the implementation fulfills all current requirements for this issue.
