import numpy as np


EARTH_RADIUS = 6371 * 10^3 # m


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
