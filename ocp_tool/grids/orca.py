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
            return nc.variables[f'gphi{subgrid}'][0, :, :].T.data

    def cell_longitudes(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            return nc.variables[f'glam{subgrid}'][0, :, :].T.data

    def cell_areas(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            return \
                nc.variables[f'e1{subgrid}'][0, :, :].T.data \
                * nc.variables[f'e2{subgrid}'][0, :, :].T.data

    def cell_masks(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.file, mode='r') as nc:
            if subgrid == 't':
                return nc.variables['top_level'][0, :, :].T.data
            else:
                raise NotImplementedError(
                    'Masks only implemented for T-grid so far'
                )

    def cell_corners(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')

        with Dataset(self.file, mode='r') as nc:
            '''
            Corner layout:
            i
            ^  2 ------- 1
            |  |         |
            |  |         |
            |  |         |
            |  3 --------4
            |
            +------------> j
            '''
            if subgrid == 't':
                lats = nc.variables['gphif'][0, :, :].T.data
                lons = nc.variables['glamf'][0, :, :].T.data
            elif subgrid == 'u':
                lats = nc.variables['gphiv'][0, :, :].T.data
                lons = nc.variables['glamv'][0, :, :].T.data
            elif subgrid == 'v':
                lats = nc.variables['gphiu'][0, :, :].T.data
                lons = nc.variables['glamu'][0, :, :].T.data
            else:
                raise NotImplementedError(
                     'ORCA corners implemented for t/u/v-grids only'
                )

        assert lats.shape == lons.shape

        corners = np.zeros((2, 4, *lats.shape))
        corners[:, 0, :, :] = [lats[:, :], lons[:, :]]
        corners[:, 1, 1:, :] = [lats[:-1, :], lons[:-1, :]]
        corners[:, 2, 1:, 1:] = [lats[:-1, :-1], lons[:-1, :-1]]
        corners[:, 3, :, 1:] = [lats[:, :-1], lons[:, :-1]]

        corners[:, 1, 0, :] = corners[:, 1, -1, :]
        corners[:, 2, 0, :] = corners[:, 2, -1, :]
        corners[:, 2, :, 0] = corners[:, 2, :, -1]
        corners[:, 3, :, 0] = corners[:, 3, :, -1]
        return corners
