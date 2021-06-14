import numpy as np

from .gaussian import EARTH_RADIUS


def _equidistant(start, end, N):
    '''Divides the interval [start, end]  into N subintervals of equal size and
    returns the midpoints of the subintervals'''
    return np.linspace(start, end, 2*N+1)[1::2]


def _interval_bounds(start, centers, end, *, loc, wrap=False):
    center_midpoints = 0.5*(centers[:-1]+centers[1:])
    if loc in ('l', 'left', 'lower'):
        return np.array((start if not wrap else end, *center_midpoints))
    elif loc in ('r', 'right', 'u', 'upper'):
        return np.array((*center_midpoints, end if not wrap else start))
    raise ValueError(f"Invalid value for 'loc' argument: {loc}")


def _row_distribute(row, nrows):
    return np.tile(row, (nrows, 1))


def _col_distribute(col, ncols):
    return np.tile(col, (ncols, 1)).T


class RegularLatLonGrid:

    def __init__(self, nlat, nlon, lats_start=-90):
        self._OP = lats_start
        self.lats = _equidistant(self._OP, -self._OP, nlat)
        self.lons = _equidistant(0, 360, nlon)

    @property
    def nlat(self):
        return len(self.lats)

    @property
    def nlon(self):
        return len(self.lons)

    def cell_latitudes(self):
        return _col_distribute(self.lats, len(self.lons))

    def cell_longitudes(self):
        return _row_distribute(self.lons, len(self.lats))

    def _cell_corner_latitudes(self):
        upper_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='u')
        lower_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='l')
        return np.array(
            [
                _col_distribute(upper_lats, len(self.lons)),  # corner 1
                _col_distribute(upper_lats, len(self.lons)),  # corner 2
                _col_distribute(lower_lats, len(self.lons)),  # corner 3
                _col_distribute(lower_lats, len(self.lons)),  # corner 4
            ]
        )

    def _cell_corner_longitudes(self):
        left_lons = _interval_bounds(0, self.lons, 360, loc='l')
        right_lons = _interval_bounds(0, self.lons, 360, loc='r', wrap=True)
        return np.array(
            [
                _row_distribute(right_lons, len(self.lats)),  # corner 1
                _row_distribute(left_lons, len(self.lats)),  # corner 2
                _row_distribute(left_lons, len(self.lats)),  # corner 3
                _row_distribute(right_lons, len(self.lats)),  # corner 4
            ]
        )

    def cell_corners(self):
        return np.array(
            [self._cell_corner_latitudes(), self._cell_corner_longitudes()]
        )

    def cell_areas(self):
        upper_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='u')
        lower_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='l')
        return _col_distribute(
            2*np.pi*EARTH_RADIUS**2
            * np.abs(
                np.sin(np.radians(upper_lats))
                - np.sin(np.radians(lower_lats))
            )/len(self.lons),
            len(self.lons)
        )
