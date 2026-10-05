import pytest
import math
from pathlib import Path
import pytransdatro 

@pytest.fixture
def st70_pnts():
    """Dictionary of points with Stereo70 coordinates for testing.
    These coordinates were extracted from the 
    "Help_TransDatRO_code_source_EN.pdf" document found in 
    "TransDatRO_code_source_1.03" folder at link: 
    https://rompos.ro/index.php/download/category/2-software
    """
    return {
        'P1': (693771.731, 310723.518, 122.714),
        'P2': (721361.806, 641283.450, 217.451),
        'P3': (516470.189, 165265.572,  86.267), 
        'P4': (402327.815, 713143.130,  22.941), 
        'P5': (329703.378, 333185.413, 260.515),
        'P6': (249343.594, 518651.464,  89.294), 
        'P7': (528076.247, 411159.899, 494.894),
        'P8': (334634.564, 593783.040, 100), # point in the 1D grid of Bucharest area
        'P9': (340134.564, 577283.040, 100) # point next to 1D grid of Bucharest area
    }

@pytest.fixture
def etrs89_pnts():
    """Dictionary of points with ETRS89 coordinates in radians for testing.
    These coordinates were computed using "TransDatRO v4.08", the official 
    application for transformation between Stereo70 and ETRS89 coordinate 
    reference systems. Download link is available here: 
    https://rompos.ro/index.php/download/category/2-software
    The DMS result of TransDatRO v4.06 is included below:
        P1, 47°42'56.40000"N, 22°28'31.99998"E,   162.016
        P2, 47°58'33.20000"N, 26°53'26.70002"E,   250.709
        P3, 46°03'57.39999"N, 20°40'11.60000"E,   129.254
        P4, 45°05'18.20001"N, 27°42'23.99999"E,    54.842
        P5, 44°26'51.30001"N, 22°54'09.30000"E,   301.996
        P6, 43°44'37.19999"N, 25°13'48.10001"E,   128.748
        P7, 46°14'47.59999"N, 23°50'46.10001"E,   535.707   
        P8, 44°30'19.18512"N, 26°10'40.57250"E,   135.106 
        P9, 44°33'24.48394"N, 25°58'16.63147"E,   134.991
    """
    return {
        'P1': (0.832795488117441, 0.392272445562385, 162.016),     
        'P2': (0.8373372226820751, 0.4693321259276279, 250.709),  
        'P3': (0.8040024035478642, 0.3607576171325035, 129.254), 
        'P4': (0.7869408405792202, 0.4835725580374142, 54.842),   
        'P5': (0.7757566737697043, 0.39972548637904465, 301.996), 
        'P6': (0.7634710101797448, 0.44034705514033184, 128.748),  
        'P7': (0.8071546621024385, 0.4161936375474547, 535.707),
        'P8': (0.7767645292239739, 0.4568911886359511, 135.106),
        'P9': (0.7776628832542685, 0.4532844607431238, 134.991) 
    } 

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
def st70_coo_file():
    """File name for Stereo70 coordinates file used for testing
    against the TransDatRO results. The format of the file
    is compatible with TransDatRo requirements.
    The file is expected to be found in the 'tests/data' folder.
    """
    return Path(__file__).parent / "data" / "st70_to_etrs89_input.csv"    


@pytest.fixture
def etrs89_coo_file():
    """File name for ETRS89 coordinates file used for testing
    against the TransDatRO results. The format of the file
    is compatible with TransDatRo requirements.
    The file is expected to be found in the 'tests/data' folder.
    """
    return Path(__file__).parent / "data" / "st70_to_etrs89_expected.csv"

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



def test_st70_to_etrs89_2D_fromfile(st70_coo_file, etrs89_coo_file, coo_rad_tol):
    """Test for the Stereo70 to ETRS89 2D coordinate transformation 
    by using test coordinates files computed with TransDatRO.
    """
    # arrannge
    t = pytransdatro.TransRO()
    tst_coo_st70 = []
    tst_coo_etrs89 = []
    sut = []

    with open(st70_coo_file, 'r') as f_st70:
        next(f_st70) # skip first line
        for line in f_st70:
            vals = line.split(',')
            tst_coo_st70.append((float(vals[1]), float(vals[2])))

    with open(etrs89_coo_file, 'r',  encoding='ANSI') as f_etrs89:
        next(f_etrs89) #skip  first line
        for line in f_etrs89:
            vals = [v.strip() for v in line.split(',')]
            if vals[0] == 'ENDF':
                break
            else:
                tst_coo_etrs89.append((
                    pytransdatro.utils.sexa_dms_chars_to_rad(vals[1]),
                    pytransdatro.utils.sexa_dms_chars_to_rad(vals[2])
                ))
    
    if len(tst_coo_st70) != len(tst_coo_etrs89):
        raise ValueError(f'Coordinates count in Stereo70 {len(tst_coo_st70)} is different from the ones in ETRS89 {len(tst_coo_etrs89)}.')

    # act
    for coo_st70 in tst_coo_st70:
        sut.append(t.st70_to_etrs89(coo_st70[0], coo_st70[1]))

    # assert
    for i, coo_etrs89 in enumerate(tst_coo_etrs89):
        assert math.isclose(sut[i][0], coo_etrs89[0], abs_tol = coo_rad_tol)
        assert math.isclose(sut[i][1], coo_etrs89[1], abs_tol = coo_rad_tol)

POINT_IDS = [f"P{i}" for i in range(1, 10)]

@pytest.mark.parametrize("point_id", POINT_IDS)
def test_st70_to_etrs89_2D(point_id, st70_pnts, etrs89_pnts, coo_rad_tol):
    """Test for the Stereo 70 to ETRS89 coordinate transformation without
    elevation
    (N,E) -> (lat,lon)
    Input of "st70_pnts" should transform to "etrs89_pnts"
    """    
    # arrannge
    t = pytransdatro.TransRO()
    
    # act
    lat, lon = t.st70_to_etrs89(st70_pnts[point_id][0], st70_pnts[point_id][1])
    
    # assert
    assert math.isclose(lat, etrs89_pnts[point_id][0], abs_tol = coo_rad_tol)
    assert math.isclose(lon, etrs89_pnts[point_id][1], abs_tol = coo_rad_tol)

@pytest.mark.parametrize("point_id", POINT_IDS)
def test_st70_to_etrs89_3D(point_id, st70_pnts, etrs89_pnts, coo_rad_tol, coo_elev_tol):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with elevation
    (N,E,H) -> (lat,lon,h)
    Input of "st70_pnts" should transform to "etrs89_pnts" (elevation included)
    """    
    # arrannge
    t = pytransdatro.TransRO()
    
    # act
    lat, lon, h = t.st70_to_etrs89(st70_pnts[point_id][0], st70_pnts[point_id][1], st70_pnts[point_id][2])
    
    # assert
    assert math.isclose(lat, etrs89_pnts[point_id][0], abs_tol = coo_rad_tol)
    assert math.isclose(lon, etrs89_pnts[point_id][1], abs_tol = coo_rad_tol)
    assert math.isclose(h, etrs89_pnts[point_id][2], abs_tol = coo_elev_tol)

@pytest.mark.parametrize("point_id", POINT_IDS)
def test_etrs89_to_st70_2D(point_id, st70_pnts, etrs89_pnts, coo_plan_tol):
    """Test for the ETRS89 to Stereo 70 coordinate transformation without
    elevation
    (lat,lon) -> (N,E)
    Input of "etrs89_pnts" should transform to "st70_pnts"
    """    
    # arrannge
    t = pytransdatro.TransRO()
    
    # act
    n, e = t.etrs89_to_st70(etrs89_pnts[point_id][0], etrs89_pnts[point_id][1])
    
    # assert
    assert math.isclose(n, st70_pnts[point_id][0], abs_tol = coo_plan_tol)
    assert math.isclose(e, st70_pnts[point_id][1], abs_tol = coo_plan_tol)

@pytest.mark.parametrize("point_id", POINT_IDS)
def test_etrs89_to_st70_3D(point_id, st70_pnts, etrs89_pnts, coo_plan_tol, coo_elev_tol):
    """Test for the ETRS89 to Stereo 70 coordinate transformation with elevation
    (lat,lon, h) -> (N,E,Z)
    Input of "etrs89_pnts" should transform to "st70_pnts" (elevation included)
    """    
    # arrannge
    t = pytransdatro.TransRO()
    
    # act
    n, e, z = t.etrs89_to_st70(etrs89_pnts[point_id][0], etrs89_pnts[point_id][1], etrs89_pnts[point_id][2])
    
    # assert
    assert math.isclose(n, st70_pnts[point_id][0], abs_tol = coo_plan_tol)
    assert math.isclose(e, st70_pnts[point_id][1], abs_tol = coo_plan_tol)  
    assert math.isclose(z, st70_pnts[point_id][2], abs_tol = coo_elev_tol)  

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
    # arrannge
    t = pytransdatro.TransRO() 
    
    # act
    with pytest.raises(pytransdatro.exceptions.OutOfGridErr):
        t.st70_to_etrs89(n, e)

@pytest.mark.parametrize("n", [224634.565, 235634.564,763634.564, 774634.563])
@pytest.mark.parametrize("e", [120783.041, 131783.040, 868783.040, 879783.039])
def test_st70_to_etrs89_no_grid_inner_no_data(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates where grid has no data available. 
    Input locations: 
        - coordinates 1 mm away (towards center of the grid) from the grid 
        corner ± 1 x grid step
        - grid corner ± 2 x grid step
        Should raise NoDataGridErr exception.
    """   
    # arrannge
    t = pytransdatro.TransRO() 
    
    # act
    with pytest.raises(pytransdatro.exceptions.NoDataGridErr):
        t.st70_to_etrs89(n, e)

@pytest.mark.parametrize(
    'n, e', [(549364, 147553), (736151, 477652),
             (549514, 763844), (252307, 720426)]
)
def test_st70_to_etrs89_just_out_ro_border_no_data(n, e):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with 
    coordinates out of the grid. Grid extents provided as input. 
    Should raise OutOfGridErr exception
    """   
    # arrannge
    t = pytransdatro.TransRO() 
    
    # act
    with pytest.raises(pytransdatro.exceptions.NoDataGridErr):
        t.st70_to_etrs89(n, e)