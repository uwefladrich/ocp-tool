import numpy as np

from .gaussian import EARTH_RADIUS


class RegularLatLonGrid:

    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny

    def cell_latitudes(self):
        return np.tile(
            np.linspace(90.0, -90.0, 2*self.ny, endpoint=False)[1::2],
            (self.nx, 1)
        ).T

    def cell_longitudes(self):
        return np.tile(
            np.linspace(0, 360, self.nx+1)[:-1],
            (self.ny, 1)
        )

    def cell_corners(self):
        border_lats = np.linspace(90.0, -90.0, self.ny+1)
        c_lat = np.empty((4, self.ny, self.nx))
        c_lat[0, :, :] = np.tile(border_lats[:-1], (self.nx, 1)).T  # north
        c_lat[1, :, :] = np.tile(border_lats[:-1], (self.nx, 1)).T  # north
        c_lat[2, :, :] = np.tile(border_lats[1:], (self.nx, 1)).T  # south
        c_lat[3, :, :] = np.tile(border_lats[1:], (self.nx, 1)).T  # south

        border_lons = np.linspace(0, 360.0, 2*self.nx+1)[1::2]
        border_lons = np.insert(border_lons, 0, -border_lons[0])
        c_lon = np.empty((4, self.ny, self.nx))
        c_lon[0, :, :] = np.tile(border_lons[1:], (self.ny, 1))  # east
        c_lon[1, :, :] = np.tile(border_lons[:-1], (self.ny, 1))  # west
        c_lon[2, :, :] = np.tile(border_lons[:-1], (self.ny, 1))  # east
        c_lon[3, :, :] = np.tile(border_lons[1:], (self.ny, 1))  # west

        return np.array([c_lat, c_lon])

    def cell_areas(self):
        sine = np.sin(
            np.radians(
                np.linspace(90.0, -90.0, 2*self.ny, endpoint=False)[1::2]
            )
        )
        return np.tile(
            2*np.pi*EARTH_RADIUS**2*(sine[:-1]-sine[1:]) / self.nx,
            (self.nx, 1)
        )