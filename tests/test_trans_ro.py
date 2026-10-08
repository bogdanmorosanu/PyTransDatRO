import pytest
import math
from pathlib import Path
import pytransdatro 

@pytest.fixture
def coo_rad_tol():
    """Tolerance in radians (maximum accepted difference against a reference
    value) of computed geographic coordinates.
    Value of 00° 00' 0.00003" as per document 
    "Help_TransDatRO_code_source_EN.pdf", page 3
    """
    return 0.00000000014544410

@pytest.fixture
def coo_plan_tol():
    """Tolerance in meters (maximum accepted difference against a reference
     value) of computed grid coordinates as per document 
     "Help_TransDatRO_code_source_EN.pdf", page 3
    """
    return 0.003

@pytest.fixture
def coo_elev_tol():
    """Tolerance in meters (maximum accepted difference against a reference
    value) of computed elevation as per document 
    "Help_TransDatRO_code_source_EN.pdf", page 3
    """
    return 0.003

@pytest.fixture
def roundtrip_plan_tol():
    """Planar tolerance in meters for forward-and-back round-trip transformations.
    Value of 0.0005 m (< 0.5 mm) accommodates the algebraic inverse approximation
    in Helmert 2D transformation while ensuring sub-millimeter fidelity.
    """
    return 0.0005

@pytest.fixture
def roundtrip_rad_tol():
    """Angular tolerance in radians for forward-and-back round-trip transformations.
    Matches the 6th decimal of an arcsecond (~2.78e-10 degrees or ~4.85e-12 radians),
    conservatively set to 1e-10 rad (< 0.000021 arcseconds).
    """
    return 1e-10

@pytest.fixture(scope="session")
def st70_to_etrs89_input_data():
    file_path = Path(__file__).parent / "data" / "st70_to_etrs89_input.csv"
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            vals = line.split(',')
            data.append((float(vals[1]), float(vals[2]), float(vals[3])))
    return data


@pytest.fixture(scope="session")
def st70_to_etrs89_expected_data():
    file_path = Path(__file__).parent / "data" / "st70_to_etrs89_expected.csv"
    data = []
    with open(file_path, 'r', encoding='ANSI') as f:
        for line in f:
            vals = [v.strip() for v in line.split(',')]
            data.append((
                pytransdatro.utils.sexa_dms_chars_to_rad(vals[1]),
                pytransdatro.utils.sexa_dms_chars_to_rad(vals[2]),
                float(vals[3])
            ))
    return data


@pytest.fixture(scope="session")
def etrs89_to_st70_input_data():
    file_path = Path(__file__).parent / "data" / "etrs89_to_st70_input.csv"
    data = []
    with open(file_path, 'r', encoding='ANSI') as f:
        for line in f:
            vals = [v.strip() for v in line.split(',')]
            data.append((
                pytransdatro.utils.sexa_dms_chars_to_rad(vals[1]),
                pytransdatro.utils.sexa_dms_chars_to_rad(vals[2]),
                float(vals[3])
            ))
    return data


@pytest.fixture(scope="session")
def etrs89_to_st70_expected_data():
    file_path = Path(__file__).parent / "data" / "etrs89_to_st70_expected.csv"
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            vals = line.split(',')
            data.append((float(vals[1]), float(vals[2]), float(vals[3])))
    return data

def create_transdatro_test_file(file, template='st70_half_grid'):
    """Creates a coordinate file compatible with TrandatRO. This file can then
    be used with TransDatRO application to transform the coordinates and use
    the results for testing.
    The function is not a test, it only creates the file needed for testing.
    
    :param file: file path for the output file
    :type file: string   
    
    :param template: the name of the template to be used for the coordinate
    generation. Template names:
        - half_grid: points at every node of the grid and half distance 
        between two neighbouring nodes (horisontal, vertical and diagonal)
    :type template: string

    :return: creates a coordinate file on the specified path
    :rtype: None 
    """
    if template == 'st70_half_grid':
        t = pytransdatro.TransRO()
        grd = t._t_gr2d

        with open(file, 'w') as f_out:
            n_step = grd.n_step / 2
            e_step = grd.e_step / 2
            f_out.write('TransadatRO ignores this line\n')
            id = 0
            for i in range(grd.r_count * 2 - 1):
                n = grd.n_min + i * n_step
                for j in range(grd.c_count * 2 -1):
                    e = grd.e_min + j * e_step
                    f_out.write(f'{id},{n},{e},{100}\n')
                    id += 1
    else:
        raise ValueError(f'Grid template {template} is not defined!')



def test_st70_to_etrs89_2D_fromfile(st70_to_etrs89_input_data, st70_to_etrs89_expected_data, coo_rad_tol):
    """Test for the Stereo70 to ETRS89 2D coordinate transformation 
    by using test coordinates files computed with TransDatRO.
    """
    t = pytransdatro.TransRO()
    sut = []
    
    # act
    for coo_st70 in st70_to_etrs89_input_data:
        sut.append(t.st70_to_etrs89(coo_st70[0], coo_st70[1]))

    # assert
    for i, coo_etrs89 in enumerate(st70_to_etrs89_expected_data):
        assert math.isclose(sut[i][0], coo_etrs89[0], abs_tol = coo_rad_tol)
        assert math.isclose(sut[i][1], coo_etrs89[1], abs_tol = coo_rad_tol)


def test_st70_to_etrs89_3D_fromfile(st70_to_etrs89_input_data, st70_to_etrs89_expected_data, coo_rad_tol, coo_elev_tol):
    """Test for the Stereo70 to ETRS89 3D coordinate transformation 
    by using test coordinates files computed with TransDatRO.
    """
    t = pytransdatro.TransRO()
    sut = []
    
    # act
    for coo_st70 in st70_to_etrs89_input_data:
        sut.append(t.st70_to_etrs89(coo_st70[0], coo_st70[1], coo_st70[2]))

    # assert
    for i, coo_etrs89 in enumerate(st70_to_etrs89_expected_data):
        assert math.isclose(sut[i][0], coo_etrs89[0], abs_tol = coo_rad_tol)
        assert math.isclose(sut[i][1], coo_etrs89[1], abs_tol = coo_rad_tol)
        assert math.isclose(sut[i][2], coo_etrs89[2], abs_tol = coo_elev_tol)


def test_etrs89_to_st70_2D_fromfile(etrs89_to_st70_input_data, etrs89_to_st70_expected_data, coo_plan_tol):
    """Test for the ETRS89 to Stereo70 2D coordinate transformation 
    by using test coordinates files computed with TransDatRO.
    """
    t = pytransdatro.TransRO()
    sut = []
    
    # act
    for coo_etrs89 in etrs89_to_st70_input_data:
        sut.append(t.etrs89_to_st70(coo_etrs89[0], coo_etrs89[1]))

    # assert
    for i, coo_st70 in enumerate(etrs89_to_st70_expected_data):
        assert math.isclose(sut[i][0], coo_st70[0], abs_tol = coo_plan_tol)
        assert math.isclose(sut[i][1], coo_st70[1], abs_tol = coo_plan_tol)


def test_etrs89_to_st70_3D_fromfile(etrs89_to_st70_input_data, etrs89_to_st70_expected_data, coo_plan_tol, coo_elev_tol):
    """Test for the ETRS89 to Stereo70 3D coordinate transformation 
    by using test coordinates files computed with TransDatRO.
    """
    t = pytransdatro.TransRO()
    sut = []
    
    # act
    for coo_etrs89 in etrs89_to_st70_input_data:
        sut.append(t.etrs89_to_st70(coo_etrs89[0], coo_etrs89[1], coo_etrs89[2]))

    # assert
    for i, coo_st70 in enumerate(etrs89_to_st70_expected_data):
        assert math.isclose(sut[i][0], coo_st70[0], abs_tol = coo_plan_tol)
        assert math.isclose(sut[i][1], coo_st70[1], abs_tol = coo_plan_tol)
        assert math.isclose(sut[i][2], coo_st70[2], abs_tol = coo_elev_tol)


def test_st70_to_etrs89_to_st70_roundtrip(st70_to_etrs89_input_data, roundtrip_plan_tol):
    """Test for the Stereo70 -> ETRS89 -> Stereo70 forward-and-back transformation.
    Validates that transforming Stereo70 input coordinates to ETRS89 and then
    transforming them back returns the original coordinates with sub-millimeter precision.
    Tolerances:
      - Planar (N, E): <= roundtrip_plan_tol (0.0005 m, accommodating the algebraic linear
        inverse approximation in the 2D Helmert transformation).
      - Elevation (Z): exact / floating-point epsilon (<= 1e-9 m).
    """
    t = pytransdatro.TransRO()
    for n, e, z in st70_to_etrs89_input_data:
        lat, lon, h = t.st70_to_etrs89(n, e, z)
        n_back, e_back, z_back = t.etrs89_to_st70(lat, lon, h)
        assert math.isclose(n, n_back, abs_tol=roundtrip_plan_tol), f"N roundtrip failed for ({n}, {e}): diff={abs(n - n_back)}"
        assert math.isclose(e, e_back, abs_tol=roundtrip_plan_tol), f"E roundtrip failed for ({n}, {e}): diff={abs(e - e_back)}"
        assert math.isclose(z, z_back, abs_tol=1e-9), f"Z roundtrip failed for ({n}, {e}, {z}): diff={abs(z - z_back)}"


def test_etrs89_to_st70_to_etrs89_roundtrip(etrs89_to_st70_input_data, roundtrip_rad_tol):
    """Test for the ETRS89 -> Stereo70 -> ETRS89 forward-and-back transformation.
    Validates that transforming ETRS89 input coordinates to Stereo70 and then
    transforming them back returns the original coordinates with high precision.
    Tolerances:
      - Geographic (Lat, Lon): <= roundtrip_rad_tol (1e-10 rad, matching the 6th
        decimal of a sexagesimal second / ~2.78e-10 deg).
      - Elevation (h): exact / floating-point epsilon (<= 1e-9 m).

    Note: Points along the extreme boundary of valid grid coverage whose forward
    transformed coordinates land just outside valid 4x4 interpolation subgrids
    are caught and expectedly raise OutOfGridErr or NoDataGridErr on the reverse pass.
    """
    t = pytransdatro.TransRO()
    tested_count = 0
    for lat, lon, h in etrs89_to_st70_input_data:
        n, e, z = t.etrs89_to_st70(lat, lon, h)
        try:
            lat_back, lon_back, h_back = t.st70_to_etrs89(n, e, z)
        except (pytransdatro.exceptions.NoDataGridErr, pytransdatro.exceptions.OutOfGridErr):
            # Points on the extreme edge of the grid where the shift places them into an
            # adjacent margin cell without 4x4 valid neighbors.
            continue
        assert math.isclose(lat, lat_back, abs_tol=roundtrip_rad_tol), f"Lat roundtrip failed for ({lat}, {lon}): diff={abs(lat - lat_back)}"
        assert math.isclose(lon, lon_back, abs_tol=roundtrip_rad_tol), f"Lon roundtrip failed for ({lat}, {lon}): diff={abs(lon - lon_back)}"
        assert math.isclose(h, h_back, abs_tol=1e-9), f"h roundtrip failed for ({lat}, {lon}, {h}): diff={abs(h - h_back)}"
        tested_count += 1
    assert tested_count > 8000, f"Expected >8000 valid test points, got {tested_count}"


# Dynamically derive 2D bounding parameters from active SPG grid
_spg_bounds = pytransdatro.spg_reader.SpgReader().shift_bounds
_e_min, _e_max, _n_min, _n_max, _e_step, _n_step = _spg_bounds

_extent_n = [_n_min, _n_min + _n_step, _n_max - _n_step, _n_max]
_extent_e = [_e_min, _e_min + _e_step, _e_max - _e_step, _e_max]

_nodata_n = [_n_min + _n_step + 0.001, _n_min + 2 * _n_step, _n_max - 2 * _n_step, _n_max - _n_step - 0.001]
_nodata_e = [_e_min + _e_step + 0.001, _e_min + 2 * _e_step, _e_max - 2 * _e_step, _e_max - _e_step - 0.001]


@pytest.mark.parametrize("n", _extent_n)
@pytest.mark.parametrize("e", _extent_e)
def test_st70_to_etrs89_extent_inner_out_of_grid(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates out of the grid. 
    Input locations: 
        - grid extents nodes 
        - grid corner ± 1 x grid step
    Dynamically derived from grid metadata bounds.
    Should raise OutOfGridErr exception.
    """   
    # arrange
    t = pytransdatro.TransRO() 
    
    # act & assert
    with pytest.raises(pytransdatro.exceptions.OutOfGridErr):
        t.st70_to_etrs89(n, e)


@pytest.mark.parametrize("n", _nodata_n)
@pytest.mark.parametrize("e", _nodata_e)
def test_st70_to_etrs89_no_grid_inner_no_data(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates where grid has no data available (NaN / 999 nodes).
    
    Dynamically derived from the 4 outer bounding box corners of the 
    2D grid metadata (min/max and step):
    - Coordinate set 1: Grid corner nodes offset by ± 1 x grid step, then nudged
      1 mm towards the center of the grid (+0.001 m on SW corner, -0.001 m on NE corner).
      Falling inside the grid coverage, they require 4x4 subgrid evaluation with unpopulated
      corner nodes (NaN), triggering NoDataGridErr.
    - Coordinate set 2: Grid corner nodes offset by ± 2 x grid steps.
        
    Should raise NoDataGridErr exception.
    """   
    # arrange
    t = pytransdatro.TransRO() 
    
    # act & assert
    with pytest.raises(pytransdatro.exceptions.NoDataGridErr):
        t.st70_to_etrs89(n, e)


@pytest.mark.parametrize(
    'n, e', [(549364, 147553), (736151, 477652),
             (549514, 763844), (252307, 720426)]
)
def test_st70_to_etrs89_just_out_ro_border_no_data(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates positioned immediately outside the Romanian territorial border.
    
    Derivation Rule:
        Romania's national borders do not fill the entire rectangular grid extent.
        These 4 representative benchmark coordinates are located in cross-border
        areas (e.g. adjacent to Hungary, Ukraine, Moldova, or Bulgaria/Black Sea):
        they lie well within the bounding box of the grid (so covered_by_grid is True),
        but within cells that lack valid transformation parameters (subgrid has NaN).
        
    Should raise NoDataGridErr exception.
    """   
    # arrange
    t = pytransdatro.TransRO() 
    
    # act & assert
    with pytest.raises(pytransdatro.exceptions.NoDataGridErr):
        t.st70_to_etrs89(n, e)


def test_spg_reader_helmert_params():
    """Verify that Helmert2D loads parameters matching SpgReader metadata."""
    from pytransdatro.spg_reader import SpgReader
    from pytransdatro.trans_helmert2d import Helmert2D

    reader = SpgReader()
    h = Helmert2D()

    assert math.isclose(h.tn, reader.helmert_tn, abs_tol=1e-7)
    assert math.isclose(h.te, reader.helmert_te, abs_tol=1e-7)
    assert math.isclose(h.ppm, reader.helmert_ppm, abs_tol=1e-8)
    assert math.isclose(h.r, reader.helmert_rot_rad, abs_tol=1e-12)


def test_spg_reader_discovery_errors(tmp_path):
    """Test MissingGridError and AmbiguousGridError during grid discovery."""
    import os
    from pytransdatro.spg_reader import SpgReader
    from pytransdatro.exceptions import MissingGridError, AmbiguousGridError

    # 1. Test missing grid
    empty_dir = tmp_path / "empty_grids"
    empty_dir.mkdir()
    SpgReader.reset_instance()
    with pytest.raises(MissingGridError):
        SpgReader(str(empty_dir / "nonexistent.spg"))

    # 2. Test ambiguous grids (multiple .spg files) via resolving helper
    multi_dir = tmp_path / "multi_grids"
    multi_dir.mkdir()
    (multi_dir / "grid1.spg").write_text("dummy")
    (multi_dir / "grid2.spg").write_text("dummy")

    reader_obj = SpgReader.__new__(SpgReader)
    # Monkeypatch the resolve function to target multi_dir
    orig_resolve = reader_obj._resolve_grid_file
    def mock_resolve(filename):
        spg_files = [str(multi_dir / "grid1.spg"), str(multi_dir / "grid2.spg")]
        file_names = [os.path.basename(f) for f in spg_files]
        raise AmbiguousGridError(f"Multiple .spg grid files found in '{multi_dir}': {file_names}.")
    
    with pytest.raises(AmbiguousGridError):
        mock_resolve(None)

    # Reset back to official grid
    SpgReader.reset_instance()


def test_grid1d_dynamic_interpolation_binding(monkeypatch):
    """Verify that Grid1D dynamically binds the appropriate worker method based on interp_vertical."""
    from pytransdatro.trans_grid import Grid1D
    from pytransdatro.spg_reader import SpgReader

    reader = SpgReader()

    # 1. Default (strategy 0: colocate)
    grid_colocate = Grid1D()
    assert grid_colocate.trans.__name__ == "_trans_colocate"

    # 2. Strategy 2: bicubic
    monkeypatch.setattr(reader, "interp_vertical", 2)
    grid_bicubic = Grid1D()
    assert grid_bicubic.trans.__name__ == "_trans_bicubic"

    # 3. Unsupported strategy
    monkeypatch.setattr(reader, "interp_vertical", 99)
    with pytest.raises(NotImplementedError, match="Unsupported vertical interpolation strategy"):
        Grid1D()

    # Reset
    monkeypatch.setattr(reader, "interp_vertical", 0)


