import numpy as np
from netCDF4 import Dataset


def _valid_subgrid(subgrid):
    return subgrid in ('t', 'u', 'v')


class ORCA:

    _orca_names = {
        (362, 292, 75): 'ORCA1L75',
    }

    def __init__(self, file):

        self.file = file

        with Dataset(file, mode='r') as nc:
            if not {'x', 'y'}.issubset(nc.dimensions):
                raise RuntimeError(
                    'Missing dimensions in NetCDF file'
                )
            if not {
                'glamt', 'glamu', 'glamv', 'glamf',
                'gphit', 'gphiu', 'gphiv', 'gphif',
                'e1t', 'e1u', 'e1v', 'e1f',
                'e2t', 'e2u', 'e2v', 'e2f',
                'top_level',
            }.issubset(nc.variables):
                raise RuntimeError(
                    'Missing variables in NetCDF file'
                )
            try:
                self.name = self._orca_names[
                    (
                        nc.dimensions['x'].size,
                        nc.dimensions['y'].size,
                        nc.dimensions['z'].size,
                    )
                ]
            except KeyError:
                raise RuntimeError(
                    'Unknown dimensions in NEMO grid file'
                )

    def cell_latitudes(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            return nc.variables[f'gphi{subgrid}'][0, ...].data

    def cell_longitudes(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            return nc.variables[f'glam{subgrid}'][0, ...].data

    def cell_areas(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            return \
                nc.variables[f'e1{subgrid}'][0, ...].data \
                * nc.variables[f'e2{subgrid}'][0, ...].data

    def cell_masks(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            tmask = np.where(
                nc.variables['top_level'][0, ...].data == 0, 1, 0
            )
            if subgrid == 't':
                return tmask
            elif subgrid == 'u':
                return tmask \
                       * tmask.take(
                            range(1, tmask.shape[1]+1), axis=1, mode='wrap'
                         )
            elif subgrid == 'v':
                return tmask \
                       * tmask.take(
                            range(1, tmask.shape[0]+1), axis=0, mode='clip'
                         )

    def cell_corners(self, subgrid='t'):
        """For the ORCA grid and staggered subgrids, see
        NEMO reference manual, section 4 'Space Domain (DOM)'

        Corner numbering used here:
        j
        ^  1 ------- 0
        |  |         |
        |  |         |
        |  |         |
        |  2 --------3
        +------------> i
        """
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')

        with Dataset(self.file, mode='r') as nc:
            if subgrid == 't':
                lats = nc.variables['gphif'][0, ...].data
                lons = nc.variables['glamf'][0, ...].data
            elif subgrid == 'u':
                lats = nc.variables['gphiv'][0, ...].data
                lons = nc.variables['glamv'][0, ...].data
            elif subgrid == 'v':
                lats = nc.variables['gphiu'][0, ...].data
                lons = nc.variables['glamu'][0, ...].data

        if lats.shape != lons.shape:
            raise ValueError(f'Incompatible lat/lon arrays in {self.file}')

        # Note that some corner lats/lons will be left undefined (set to an
        # invalid initial value), because we do not handle the north-fold
        # (v-grid) or the southern end of the grid over the Antarctic (t,
        # u-grids).
        corners = np.full((2, 4, *lats.shape), -99999.0)

        if subgrid == 't':
            corners[0, 1, :, :] = np.roll(lats, 1, axis=1)
            corners[1, 1, :, :] = np.roll(lons, 1, axis=1)

            corners[0, 2, 1:, :] = np.roll(lats[:-1, :], 1, axis=1)
            corners[1, 2, 1:, :] = np.roll(lons[:-1, :], 1, axis=1)

        elif subgrid == 'u':
            corners[0, 1, :, :] = lats
            corners[1, 1, :, :] = lons

            corners[0, 2, 1:, :] = lats[:-1, :]
            corners[1, 2, 1:, :] = lons[:-1, :]

        elif subgrid == 'v':
            corners[0, 1, :-1, :] = np.roll(lats[1:, :], 1, axis=1)
            corners[1, 1, :-1, :] = np.roll(lons[1:, :], 1, axis=1)

            corners[0, 2, :, :] = np.roll(lats, 1, axis=1)
            corners[1, 2, :, :] = np.roll(lons, 1, axis=1)

        corners[:, 0, :, :] = np.roll(corners[:, 1, :, :], -1, axis=2)
        corners[:, 3, :, :] = np.roll(corners[:, 2, :, :], -1, axis=2)

        return corners
