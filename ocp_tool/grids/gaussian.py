import numpy as np


EARTH_RADIUS = 6371e6  # m


def latband_borders(lats):
    borders = np.empty(len(lats)+1)
    borders[0] = 90  # north pole
    borders[-1] = -90  # south pole
    borders[1:-1] = (lats[:-1]+lats[1:])/2
    return borders


def latband_areas(lats):
    sine = np.sin(
        np.radians(
            latband_borders(lats)
        )
    )
    return 2*np.pi*EARTH_RADIUS**2*(sine[:-1]-sine[1:])


def longitudes(nlons):
    return np.linspace(0, 360, nlons+1)[:-1]


def longitude_borders(nlons):
    borders = np.linspace(0, 360, 2*nlons+1)[:-1][1::2]
    borders = np.insert(borders, 0, -borders[0])
    return borders


class GaussianGrid:

    def __init__(self, lats):
        self.lats = np.array(lats)
        self.nlats = len(self.lats)
        self.nlons = 2*self.nlats

    def cell_latitudes(self):
        return np.tile(self.lats, (self.nlons, 1)).T

    def cell_longitudes(self):
        return np.tile(longitudes(self.nlons), (self.nlats, 1))

    def cell_corners(self):
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
        border_lats = latband_borders(self.lats)
        c_lat = np.empty((4, self.nlats, self.nlons))
        c_lat[0, :, :] = np.tile(border_lats[:-1], (self.nlons, 1)).T  # north
        c_lat[1, :, :] = np.tile(border_lats[:-1], (self.nlons, 1)).T  # north
        c_lat[2, :, :] = np.tile(border_lats[1:], (self.nlons, 1)).T  # south
        c_lat[3, :, :] = np.tile(border_lats[1:], (self.nlons, 1)).T  # south

        border_lons = longitude_borders(self.nlons)
        c_lon = np.empty((4, self.nlats, self.nlons))
        c_lon[0, :, :] = np.tile(border_lons[1:], (self.nlats, 1))  # east
        c_lon[1, :, :] = np.tile(border_lons[:-1], (self.nlats, 1))  # west
        c_lon[2, :, :] = np.tile(border_lons[:-1], (self.nlats, 1))  # east
        c_lon[3, :, :] = np.tile(border_lons[1:], (self.nlats, 1))  # west

        return np.array([c_lat, c_lon])

    def cell_areas(self):
        return np.tile(
            latband_areas(self.lats)/self.nlons,
            (self.nlons, 1)
        )


class ReducedGaussianGrid:

    def __init__(self, nlats, lats, nlons):
        self.nlats = nlats # int
        self.lats = lats # float[nlats]
        self.nlons = nlons # int[nlats]
        assert nlats == len(lats) == len(nlons)
        self.dlons = (360/nl for nl in nlons)
        self.size = sum(self.nlons)

    def cell_latitudes(self): # float[self.size]
        latitudes = np.array(
            [
                lat for sub in (
                    ((l,)*n for l,n in zip(self.lats, self.nlons))
                ) for lat in sub
            ],
            dtype=np.float64
        )
        assert len(latitudes) == self.size
        return latitudes

    def cell_longitudes(self): # float[self.size]
        longitudes = np.array(
            [
                360*i/nlon for nlon in self.nlons for i in range(nlon)
            ],
            dtype=np.float64
        )
        assert len(longitudes) == self.size
        return longitudes

    def cell_corners(self): # float[self.size, 4]
        corners = np.empty([2, 4, self.size])
        cell_idx = 0
        for i, (lat,nlon) in enumerate(zip(self.lats, self.nlons)):
            lat_n = (self.lats[i-1]+lat)/2 if i>0 else (lat+90)/2
            assert -90 < lat_n < 90
            lat_s = (self.lats[i+1]+lat)/2 if i<self.nlats-1 else (lat-90)/2
            assert -90 < lat_s < 90
            dlon = 180/nlon
            for lon in (i*360/nlon for i in range(nlon)):
                '''
                Corner layout:
                y
                ^  2 ---------- 1
                |  |            |
                |  |            |
                |  |            |
                |  3 -----------4
                |
                +---------------> x
                '''
                lon_e = lon+dlon if lon+dlon<=180 else lon+dlon-360
                assert -180 < lon_e <= 180
                lon_w = lon-dlon if lon-dlon<=180 else lon-dlon-360
                assert -180 < lon_w <= 180
                corners[:, 0, cell_idx] = [lat_n, lon_e]
                corners[:, 1, cell_idx] = [lat_n, lon_w]
                corners[:, 2, cell_idx] = [lat_s, lon_w]
                corners[:, 3, cell_idx] = [lat_s, lon_e]
                cell_idx += 1
        assert cell_idx == self.size
        return corners

    def cell_areas(self): # float[self.size]
        areas = np.empty([0])
        for i, (lat,nlon) in enumerate(zip(self.lats, self.nlons)):
            dlat_n = (self.lats[i-1]-lat)/2 if i>0 else 90-lat
            dlat_s = (lat-self.lats[i+1])/2 if i<self.nlats-1 else lat+90
            dlon = 180/nlon
            dx = dlon * np.pi/180 * EARTH_RADIUS * np.cos( np.pi/180 * lat)
            dy = (dlat_n+dlat_s) * np.pi/180 *EARTH_RADIUS
            areas = np.append(areas, [dx*dy]*nlon)
        assert len(areas) == self.size
        return areas
