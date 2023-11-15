"""Defines expections that are thrown be pyTransDat
"""
class OutOfGridErr(Exception):
    """Signifies that the input location (N,E) 
    is not cover by a grid / is out of a grid's min/max
    """
    def __init__(self, n, e, grid):
        """Constructor method

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid        
        """
        message = (
            f'({n}, {e}) location is out of the grid area defined by file {grid.source}. '
            f'Location must be inside the grid defined by '
            f'Min N = {grid.n_min}, Min E = {grid.e_min}, '
            f'Max N = {grid.n_max}, Max E  = {grid.e_max}.'
        )
        super().__init__(message)

class NoDataGridErr(Exception):
    """Signifies that the input location (N,E) 
    grid shif values required for interpolation have No Data value(s),
    therefore computation is not possible
    """
    def __init__(self, n, e, grid):
        """Constructor method

        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid
        """
        message = (
            f'(The subgrid of size {grid.sg_size}x{grid.sg_size} nodes required '
            f'for the computation of interpolated value at {n}, {e} location '
            f'has No Data value(s).'
        )
        super().__init__(message)  

class OutOfRangeIndexGridErr(Exception):
    """Signifies that the input index is out of the range of available
    grid indexes
    """

    def __init__(self, idx, grid):
        """Constructor method

        :param idx: index
        :type idx: int

        :param grid: Grid against which index coordinates are queried
        :type grid: Grid
        """        
        message = (
            f'(Grid {grid.file_name} has index values ranging from 0 to '
            f'{grid.c_count * grid.r_count - 1}. The requested index ({idx}) '
            f'is out of range.'
        )
        super().__init__(message)  