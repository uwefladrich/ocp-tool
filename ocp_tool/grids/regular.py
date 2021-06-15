import numpy as np

from .earth import RADIUS as EARTH_RADIUS


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


class LatLonGrid:

    def __init__(self, lats, lons, lats_start=-90):
        if not (lats_start < lats[0] <= lats[-1] < -lats_start):
            raise ValueError(
                'Inconsistent latitude/lats_start values'
            )
        self._OP = lats_start
        self.lats = lats
        self.lons = lons

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
                _col_distribute(upper_lats, len(self.lons)),  # 1 ---- 0
                _col_distribute(upper_lats, len(self.lons)),  # |      |
                _col_distribute(lower_lats, len(self.lons)),  # |      |
                _col_distribute(lower_lats, len(self.lons)),  # 2 ---- 3
            ]
        )

    def _cell_corner_longitudes(self):
        left_lons = _interval_bounds(0, self.lons, 360, loc='l')
        right_lons = _interval_bounds(0, self.lons, 360, loc='r', wrap=True)
        return np.array(
            [
                _row_distribute(right_lons, len(self.lats)),  # 1 ---- 0
                _row_distribute(left_lons, len(self.lats)),   # |      |
                _row_distribute(left_lons, len(self.lats)),   # |      |
                _row_distribute(right_lons, len(self.lats)),  # 2 ---- 3
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


class RegularLatLonGrid(LatLonGrid):
    def __init__(self, nlat, nlon, first_lat=-90):
        super().__init__(
            lats=_equidistant(first_lat, -first_lat, nlat),
            lons=_equidistant(0, 360, nlon),
            first_lat=first_lat
        )


class FullGaussianGrid(LatLonGrid):
    def __init__(self, lats, first_lat=-90):
        super().__init__(
            lats=lats,
            lons=_equidistant(0, 360, 2*len(lats)),
            first_lat=first_lat
        )
