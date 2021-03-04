import pytest
import math
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
        'P7': (528076.247, 411159.899, 494.894)
    }

@pytest.fixture
def etrs89_pnts():
    """Dictionary of points with ETRS89 coordinates in radians for testing.
    These coordinates were computed using "TransDatRO v4.06", the official 
    application for transformation between Stereo70 and ETRS89 coordinate 
    reference systems. Download link is available here: 
    https://rompos.ro/index.php/download/category/2-software
    The DMS result of TransDatRO v4.06 is included below:
        P1, 47°42'56.40000"N, 22°28'31.99998"E,   162.087
        P2, 47°58'33.20000"N, 26°53'26.70002"E,   250.711
        P3, 46°03'57.39999"N, 20°40'11.60000"E,   129.232
        P4, 45°05'18.20001"N, 27°42'23.99999"E,    54.828
        P5, 44°26'51.30001"N, 22°54'09.30000"E,   301.963
        P6, 43°44'37.19999"N, 25°13'48.10001"E,   128.724
        P7, 46°14'47.59999"N, 23°50'46.10001"E,   535.704    
    """
    return {
        'P1': (0.832795488117441,  0.392272445562385,   162.087),     
        'P2': (0.8373372226820751, 0.4693321259276279,  250.711),  
        'P3': (0.8040024035478642, 0.3607576171325035,  129.232), 
        'P4': (0.7869408405792202, 0.4835725580374142,   54.828),   
        'P5': (0.7757566737697043, 0.39972548637904465, 301.963), 
        'P6': (0.7634710101797448, 0.44034705514033184, 128.724),  
        'P7': (0.8071546621024385, 0.4161936375474547,  535.704) 
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


def test_st70_to_etrs89_2D(st70_pnts, etrs89_pnts, coo_rad_tol):
    """Test for the Stereo 70 to ETRS89 coordinate transformation without
    elevation
    (N,E) -> (lat,lon)
    Input of "st70_pnts" should transform to "etrs89_pnts"
    """    
    # arrannge
    t = pytransdatro.TransRO()
    sut = {}
    
    # act
    for k, v in st70_pnts.items():
        sut[k] = t.st70_to_etrs89(v[0], v[1])
    
    # assert
    for k in sut:
        assert math.isclose(sut[k][0], etrs89_pnts[k][0], abs_tol = coo_rad_tol)
        assert math.isclose(sut[k][1], etrs89_pnts[k][1], abs_tol = coo_rad_tol)

def test_st70_to_etrs89_3D(st70_pnts, etrs89_pnts, coo_rad_tol, coo_elev_tol):
    """Test for the Stereo 70 to ETRS89 coordinate transformation with elevation
    (N,E,H) -> (lat,lon,h)
    Input of "st70_pnts" should transform to "etrs89_pnts" (elevation included)
    """    
    # arrannge
    t = pytransdatro.TransRO()
    sut = {}
    
    # act
    for k, v in st70_pnts.items():
        sut[k] = t.st70_to_etrs89(v[0], v[1], v[2])
    
    # assert
    for k in sut:
        assert math.isclose(sut[k][0], etrs89_pnts[k][0], abs_tol = coo_rad_tol)
        assert math.isclose(sut[k][1], etrs89_pnts[k][1], abs_tol = coo_rad_tol)
        assert math.isclose(sut[k][2], etrs89_pnts[k][2], abs_tol = coo_elev_tol)

def test_etrs89_to_st70_2D(st70_pnts, etrs89_pnts, coo_plan_tol):
    """Test for the ETRS89 to Stereo 70 coordinate transformation without
    elevation
    (lat,lon) -> (N,E)
    Input of "etrs89_pnts" should transform to "st70_pnts"
    """    
    # arrannge
    t = pytransdatro.TransRO()
    sut = {}
    
    # act
    for k, v in etrs89_pnts.items():
        sut[k] = t.etrs89_to_st70(v[0], v[1])
    
    # assert
    for k in sut:
        assert math.isclose(sut[k][0], st70_pnts[k][0], abs_tol = coo_plan_tol)
        assert math.isclose(sut[k][1], st70_pnts[k][1], abs_tol = coo_plan_tol)

def test_etrs89_to_st70_3D(st70_pnts, etrs89_pnts, coo_plan_tol, coo_elev_tol):
    """Test for the ETRS89 to Stereo 70 coordinate transformation with elevation
    (lat,lon, h) -> (N,E,Z)
    Input of "etrs89_pnts" should transform to "st70_pnts" (elevation included)
    """    
    # arrannge
    t = pytransdatro.TransRO()
    sut = {}
    
    # act
    for k, v in etrs89_pnts.items():
        sut[k] = t.etrs89_to_st70(v[0], v[1], v[2])
    
    # assert
    for k in sut:
        assert math.isclose(sut[k][0], st70_pnts[k][0], abs_tol = coo_plan_tol)
        assert math.isclose(sut[k][1], st70_pnts[k][1], abs_tol = coo_plan_tol)  
        assert math.isclose(sut[k][2], st70_pnts[k][2], abs_tol = coo_plan_tol)  

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