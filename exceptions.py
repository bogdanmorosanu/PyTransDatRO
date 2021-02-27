    """Defines expections that are thrown be pyTransDat
    """
class OutOfGridErr(Exception):
    """Signifies that the input location (N,E) 
    is not cover by a grid / is out of a grid's min/maxt
    """
    def __init__(self, n, e, grid):
        """
        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid
        """
        message = (
            f'({n}, {e}) location is out of the grid area.'
            f'Location must be inside the grid defined by '
            f'Min N = {n_min}, Min E = {e_min}, '
            f'Max N = {n_max}, Max E  = {e_max}.'
        )
        super().__init__(message)

class NoDataGridErr(Exception):
    """Signifies that the input location (N,E) 
    grid shif values required for interpolation have No Data value(s),
    therefore computation is not possible
    """
    def __init__(self, n, e, grid):
        """
        :param n: northing
        :type n: float

        :param e: easting
        :type e: float 

        :param grid: Grid against which values (n, e) have been tested
        :type grid: Grid
        """
        message = (
            f'(The subgrid of size {grid.sg_size}x{grid.sg_size} required '
            f'for the computation of interpolated value at {n}, {e} location '
            f'has No Data value(s).'
        )        