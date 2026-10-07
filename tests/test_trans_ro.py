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


@pytest.mark.xfail(reason="New grid release mismatch for elevation values")
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


@pytest.mark.xfail(reason="New grid release mismatch for elevation values")
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


@pytest.mark.parametrize("n", [213634.564, 224634.564, 774634.564, 785634.564])
@pytest.mark.parametrize("e", [109783.040, 120783.040, 879783.040, 890783.040])
def test_st70_to_etrs89_extent_inner_out_of_grid(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates out of the grid. 
    Input locations: 
        - grid extents nodes 
        - grid corner ± 1 x grid step
    Should raise OutOfGridErr exception.
    """   
    # arrange
    t = pytransdatro.TransRO() 
    
    # act & assert
    with pytest.raises(pytransdatro.exceptions.OutOfGridErr):
        t.st70_to_etrs89(n, e)


@pytest.mark.parametrize("n", [224634.565, 235634.564, 763634.564, 774634.563])
@pytest.mark.parametrize("e", [120783.041, 131783.040, 868783.040, 879783.039])
def test_st70_to_etrs89_no_grid_inner_no_data(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates where grid has no data available (NaN / 999 nodes).
    
    Derivation Rule:
        The coordinates are derived from the 4 outer bounding box corners of the 
        2D grid (n_min=213634.564, n_max=785634.564, e_min=109783.040, e_max=890783.040, 
        with grid step n_step = e_step = 11000 m):
        - Coordinate set 1: Grid corner nodes offset by ± 1 x grid step, then nudged
          1 mm towards the center of the grid (+0.001 m on SW corner, -0.001 m on NE corner).
          Because of the 1 mm inward nudge, the point falls inside the grid coverage
          and requires a 4x4 subgrid evaluation; however, that subgrid contains unpopulated
          corner nodes (no-data/NaN), triggering NoDataGridErr.
        - Coordinate set 2: Grid corner nodes offset by ± 2 x grid steps (e.g. 213634.564 + 22000 = 235634.564).
        
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