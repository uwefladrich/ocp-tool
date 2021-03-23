import os

import numpy as np
from netCDF4 import Dataset as NCDataset

def get_var(nc, name, type_, dim):
    if name not in nc.variables:
        return nc.createVariable(name, type_, dim)
    else:
        return nc.variables[name]

def write_grid(name, lats, lons, corners=None, path=None, append=True):

    if lats.shape != lons.shape:
        raise ValueError('Mismatch between lat and lon dimensions')

    if lats.ndim not in (1,2):
        raise ValueError('Invalid dimensions, must be one or two dimensional')

    two_dim = lats.ndim == 2

    x_n = lats.shape[0]
    y_n = lats.shape[1] if two_dim else 1

    x_d = f'{name}_x'
    y_d = f'{name}_y'
    c_d = f'{name}_c'

    lat_v = f'{name}.lat'
    lon_v = f'{name}.lon'
    cla_v = f'{name}.cla'
    clo_v = f'{name}.clo'

    with NCDataset(
            os.path.join(path or '', 'grids.nc'),
            mode='r+' if append else 'w'
         ) as nc:

        if x_d not in nc.dimensions:
            nc.createDimension(x_d, x_n)
        if y_d not in nc.dimensions:
            nc.createDimension(y_d, y_n)

        lat_id = get_var(nc, lat_v, 'float64', (y_d, x_d))
        lat_id.units = 'degrees_north'
        lat_id.standard_name = 'Latitude'
        lat_id[:,:] = lats.T if two_dim else [lats]

        lon_id = get_var(nc, lon_v, 'float64', (y_d, x_d))
        lon_id.units = 'degrees_east'
        lon_id.standard_name = 'Longitude'
        lon_id[:,:] = lons.T if two_dim else [lons]

        if corners is not None:
            if c_d not in nc.dimensions:
                nc.createDimension(c_d, 4)

            cla_id = get_var(nc, cla_v, 'float64', (c_d, y_d, x_d))
            cla_id.units = 'degrees_north'
            cla_id.standard_name = 'Corner_latitude'

            clo_id = get_var(nc, clo_v, 'float64', (c_d, y_d, x_d))
            clo_id.units = 'degrees_east'
            clo_id.standard_name = 'Corner_longitude'

            if two_dim:
                assert corners.ndim == 4
                cla_id[:,:,:] = np.transpose(corners[0, :, :, :], axes=(0, 2, 1))
                clo_id[:,:,:] = np.transpose(corners[1, :, :, :], axes=(0, 2, 1))
            else:
                assert corners.ndim == 3
                cla_id[:,0,:] = corners[0, :, :]
                clo_id[:,0,:] = corners[1, :, :]

def write_area(name, areas, path=None, append=True):

    if areas.ndim not in (1,2):
        raise ValueError('Invalid dimensions, must be one or two dimensional')

    two_dim = areas.ndim == 2

    x_n = areas.shape[0]
    y_n = areas.shape[1] if two_dim else 1

    x_d = f'{name}_x'
    y_d = f'{name}_y'

    areas_v = f'{name}.srf'

    with NCDataset(
            os.path.join(path or '', 'areas.nc'),
            mode='r+' if append else 'w'
         ) as nc:

        if x_d not in nc.dimensions:
            nc.createDimension(x_d, x_n)
        if y_d not in nc.dimensions:
            nc.createDimension(y_d, y_n)

        areas_id = get_var(nc, areas_v, 'float64', (y_d, x_d))
        areas_id[:,:] = areas.T if two_dim else [areas]


def write_mask(name, masks, path=None, append=True):

    if masks.ndim not in (1,2):
        raise ValueError('Invalid dimensions, must be one or two dimensional')

    two_dim = masks.ndim == 2

    x_n = masks.shape[0]
    y_n = masks.shape[1] if two_dim else 1

    x_d = f'{name}_x'
    y_d = f'{name}_y'

    masks_v = f'{name}.msk'

    with NCDataset(
            os.path.join(path or '', 'masks.nc'),
            mode='r+' if append else 'w'
         ) as nc:

        if x_d not in nc.dimensions:
            nc.createDimension(x_d, x_n)
        if y_d not in nc.dimensions:
            nc.createDimension(y_d, y_n)

        masks_id = get_var(nc, masks_v, 'int32', (y_d, x_d))
        masks_id[:,:] = masks.T if two_dim else [masks]
