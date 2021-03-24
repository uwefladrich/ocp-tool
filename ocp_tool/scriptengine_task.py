import logging
import numpy as np

import ocp_tool as ocpt
from ocp_tool.grids import factory as grid_factory
from ocp_tool.oasis import write_grid, write_area, write_mask
from ocp_tool.grib import read as grib_read


try:
    from scriptengine.tasks.core import Task, timed_runner
    from scriptengine.exceptions import ScriptEngineTaskRunError
except (ModuleNotFoundError, ImportError) as e:
    logging.warning(
        'The ScriptEngine module could not be loaded, which means '
        'that the OCPTool ScriptEngine tasks will not be available. '
        f'(the error was: {e})'
    )
else:
    class OCPTool(Task):

        _required_arguments = (
            'oifs_grid_type',
            'nemo_grid_file',
        )

        def __init__(self, arguments):
            OCPTool.check_arguments(arguments)
            super().__init__(arguments)

        @timed_runner
        def run(self, context):

            self.log_info(
                'Creating OASIS files for OpenIFS and NEMO with the OCP Tool'
            )

            oifs_grid = grid_factory(
                ocpt.grids.__dict__.get(
                    self.getarg('oifs_grid_type')
                )
            )
            self.log_debug(f'OIFS grid type: {type(oifs_grid)}')

            nemo_grid = grid_factory(
                'orca',
                self.getarg('nemo_grid_file')
            )
            self.log_debug(f'NEMO grid type: {type(nemo_grid)}')

            oasis_grid_names = {
                'TQ21': 'F016',
                'TL159': 'N080',
                'TCO95': 'O096',
                'TCO159': 'O160',
                'TL255': 'N128',
                'ORCA1L75t': 'O1T0',
                'ORCA1L75u': 'O1U0',
                'ORCA1L75v': 'O1V0',
            }

            self.log_debug('Writing OASIS grids.nc file...')
            write_grid(
                name=oasis_grid_names[self.getarg('oifs_grid_type')],
                lats=oifs_grid.cell_latitudes(),
                lons=oifs_grid.cell_longitudes(),
                corners=oifs_grid.cell_corners(),
                append=False
            )
            for subgrid in ('t', 'u', 'v'):
                write_grid(
                    name=oasis_grid_names[nemo_grid.name+subgrid],
                    lats=nemo_grid.cell_latitudes(grid=subgrid),
                    lons=nemo_grid.cell_longitudes(grid=subgrid),
                    corners=nemo_grid.cell_corners(grid=subgrid)
                )
            self.log_debug('...finished writing OASIS grids.nc file')

            self.log_debug('Writing OASIS areas.nc file...')
            write_area(
                name=oasis_grid_names[self.getarg('oifs_grid_type')],
                areas=oifs_grid.cell_areas(),
                append=False
            )
            for subgrid in ('t', 'u', 'v'):
                write_area(
                    name=oasis_grid_names[nemo_grid.name+subgrid],
                    areas=nemo_grid.cell_areas(grid=subgrid)
                )
            self.log_debug('...finished writing OASIS areas.nc file')

            self.log_debug('Read and process OIFS land-sea mask')
            data = grib_read(
                self.getarg('oifs_mask_file'),
                ('lsm', 'cl')
            )
            oifs_lsm = np.where(
                np.logical_or(data['lsm']>0.5, data['cl']>0.5), 1, 0
            )
            self.log_debug('Writing OASIS masks.nc file...')
            write_mask(
                name=oasis_grid_names[self.getarg('oifs_grid_type')],
                masks=oifs_lsm,
                append=False
            )
            for subgrid in ('t', ):
                write_mask(
                    name=oasis_grid_names[nemo_grid.name+subgrid],
                    masks=nemo_grid.cell_masks(grid=subgrid)
                )
            self.log_debug('...finished writing OASIS masks.nc file')
