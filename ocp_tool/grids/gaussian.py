import numpy as np

from .earth import RADIUS as EARTH_RADIUS


def latband_borders(lats):
    '''Expects a numpy array with center latitudes and returns a numpy array
    of length len(lats)+1 with latitudes of borders between bands as mean of
    the center latitudes. North and south poles are included as first and
    last border, respectively. Everything is expected and computed in degrees
    with north pole at 90 and south pole at -90.
    '''
    borders = np.empty(len(lats)+1)
    borders[0] = 90  # north pole
    borders[-1] = -90  # south pole
    borders[1:-1] = (lats[:-1]+lats[1:])/2
    return borders


def latband_areas(lats):
    '''Expects a numpy array with center latitudes and returns a numpy array
    with areas of the latitude bands as defined by borders returned from
    latband_borders().
    '''
    sine = np.sin(
        np.radians(
            latband_borders(lats)
        )
    )
    return 2*np.pi*EARTH_RADIUS**2*(sine[:-1]-sine[1:])


def equidist_lons(nlons):
    '''Returns a numpy array with longitudes (degrees east) spaced equally
    around the globe and including the prime meridian as first element.
    '''
    return np.linspace(0, 360, nlons+1)[:-1]


def equidist_lon_borders(nlon):
    '''Returns a numpy array with nlon+1 longitude borders for longitude
    bands defined by the equidist distribution of nlons cells (see
    equidist_lons()). The first border -dlon = -360/nlon, i.e. west of the
    prime meridian.
    '''
    borders = np.linspace(0, 360, 2*nlon+1)[:-1][1::2]
    borders = np.insert(borders, 0, -borders[0])
    return borders


class ReducedGaussianGrid:
    '''An reduced Gaussian Grid, defined by a list of latitude bands and a
    number of equidistant lonitude segments for each latitude band.
    '''
    def __init__(self, lats, nlons):
        if len(lats) != len(nlons):
            raise ValueError(
                'Arrays "lats" and "nlons" do not match in size'
            )
        self.lats = np.array(lats)
        self.nlons = np.array(nlons)
        self.size = sum(self.nlons)

    def cell_latitudes(self):  # float[self.size]
        latitudes = np.empty(0)
        for lat, nlon in zip(self.lats, self.nlons):
            latitudes = np.append(latitudes, np.repeat(lat, nlon))
        assert len(latitudes) == self.size
        return latitudes

    def cell_longitudes(self):  # float[self.size]
        longitudes = np.empty(0)
        for nlon in self.nlons:
            longitudes = np.append(longitudes, equidist_lons(nlon))
        assert len(longitudes) == self.size
        return longitudes

    def cell_corners(self):  # float[self.size, 4]
        '''
        Corner layout:
        +---------------> i
        |  1 ---------- 0
        |  |            |
        |  |            |
        |  |            |
        |  2 -----------3
        v
        j
        '''
        corners = np.empty((2, 4, self.size))
        ii = 0
        brdlats = latband_borders(self.lats)
        for lat_n, lat_s, nlons in zip(brdlats[:-1], brdlats[1:], self.nlons):
            brdlons = equidist_lon_borders(nlons)
            for lon_e, lon_w in zip(brdlons[:-1], brdlons[1:]):
                corners[:, 0, ii] = [lat_n, lon_e]
                corners[:, 1, ii] = [lat_n, lon_w]
                corners[:, 2, ii] = [lat_s, lon_w]
                corners[:, 3, ii] = [lat_s, lon_e]
                ii += 1
        assert ii == self.size
        return corners

    def cell_areas(self):  # float[self.size]
        areas = np.empty([0])
        for nlon, lat_area in zip(self.nlons, latband_areas(self.lats)):
            areas = np.append(areas, np.repeat(lat_area/nlon, nlon))
        assert len(areas) == self.size
        return areas
