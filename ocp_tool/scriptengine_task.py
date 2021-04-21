import logging
import numpy as np

import ocp_tool as ocpt
import ocp_tool.grids
import ocp_tool.grib
import ocp_tool.oasis

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
            'oifs_mask_file',
        )

        def __init__(self, arguments):
            OCPTool.check_arguments(arguments)
            super().__init__(arguments)

        @timed_runner
        def run(self, context):

            oasis_grid_names = {
                'TCO159': 'ACS',  # 2nd character to be inserted later
                'TL255': 'ALS',
                'TCO95': 'ACl',
                'TL159': 'ALl',
                'TQ21': 'AQx',
                'ORCA1L75t': 'O1T0',
                'ORCA1L75u': 'O1U0',
                'ORCA1L75v': 'O1V0',
                'rnfm-atm': 'RnfA',
                'amipfr': 'AMIP',
            }

            # Construct final name for atm grid in OASIS
            def oifs_oasis_grid_name(grid_type, specifier):
                return oasis_grid_names[grid_type][0] \
                       + specifier + oasis_grid_names[grid_type][1:]

            # OpenIFS grid
            oifs_grid_type = self.getarg('oifs_grid_type', context)
            try:
                oifs_grid = ocpt.grids.factory(oifs_grid_type)
            except NotImplementedError:
                self.log_error(
                    'Invalid OIFS grid type: '
                    f'{self.getarg("oifs_grid_type", context)}'
                )
                raise ScriptEngineTaskRunError
            self.log_debug('Write OIFS grids to grids.nc')
            ocpt.oasis.write_grid(
                name=oifs_oasis_grid_name(oifs_grid_type, 'A'),
                lats=oifs_grid.cell_latitudes(),
                lons=oifs_grid.cell_longitudes(),
                corners=oifs_grid.cell_corners(),
                append=False
            )
            self.log_debug('Write OIFS areas to areas.nc')
            ocpt.oasis.write_area(
                name=oifs_oasis_grid_name(oifs_grid_type, 'A'),
                areas=oifs_grid.cell_areas(),
                append=False
            )

            oifs_mask_file = self.getarg('oifs_mask_file', context)
            try:
                oifs_masks = ocpt.grib.read(oifs_mask_file, ('lsm', 'cl'))
            except (FileNotFoundError, PermissionError):
                self.log_error(
                    f'Could not open OIFS mask file "{oifs_mask_file}"'
                )
                raise ScriptEngineTaskRunError
            oifs_lsm = np.where(
                np.logical_or(
                    oifs_masks['lsm'] > 0.5, oifs_masks['cl'] > 0.5
                ), 1, 0
            )
            self.log_debug('Write OIFS masks to masks.nc')
            ocpt.oasis.write_mask(
                name=oifs_oasis_grid_name(oifs_grid_type, 'A'),
                masks=oifs_lsm,
                append=False
            )

            # NEMO grid
            nemo_grid_file = self.getarg('nemo_grid_file', context, default=None)
            if nemo_grid_file is not None:
                try:
                    nemo_grid = ocpt.grids.factory('ORCA', nemo_grid_file)
                except (FileNotFoundError, PermissionError):
                    self.log_error(
                        f'Could not open NEMO grid file "{nemo_grid_file}"'
                    )
                    raise ScriptEngineTaskRunError
                self.log_debug('Write NEMO grids and areas to grids.nc/areas.nc')
                for subgrid in ('t', 'u', 'v'):
                    ocpt.oasis.write_grid(
                        name=oasis_grid_names[nemo_grid.name+subgrid],
                        lats=nemo_grid.cell_latitudes(subgrid=subgrid),
                        lons=nemo_grid.cell_longitudes(subgrid=subgrid),
                        corners=nemo_grid.cell_corners(subgrid=subgrid)
                    )
                    ocpt.oasis.write_area(
                        name=oasis_grid_names[nemo_grid.name+subgrid],
                        areas=nemo_grid.cell_areas(subgrid=subgrid)
                    )
                self.log_debug('Write NEMO masks and areas to masks.nc')
                for subgrid in ('t', ):
                    ocpt.oasis.write_mask(
                        name=oasis_grid_names[nemo_grid.name+subgrid],
                        masks=nemo_grid.cell_masks(subgrid=subgrid)
                    )

            # Runoff-mapper grid
            rnfm_grid = ocpt.grids.factory('F128')
            self.log_debug('Write RNFM grid to grids.nc')
            ocpt.oasis.write_grid(
                name=oasis_grid_names['rnfm-atm'],
                lats=rnfm_grid.cell_latitudes(),
                lons=rnfm_grid.cell_longitudes(),
                corners=rnfm_grid.cell_corners()
            )
            self.log_debug('Write RNFM areas to areas.nc')
            ocpt.oasis.write_area(
                name=oasis_grid_names['rnfm-atm'],
                areas=rnfm_grid.cell_areas()
            )
            self.log_debug('Write RNFM mask to masks.nc')
            if self.getarg('rnfm_mask_file', context, default=None):
                self.log_error(
                    'Reading the RNFM mask from file is not implemented yet'
                )
                raise ScriptEngineTaskRunError
            else:
                ocpt.oasis.write_mask(
                    name=oasis_grid_names['rnfm-atm'],
                    masks=np.zeros((rnfm_grid.nlons, rnfm_grid.nlats))
                )

            # AMIP Forcing-reader grid
            amipfr_grid = ocpt.grids.factory('regular_latlon', 360, 180)
            self.log_debug('Write AMIP-FR grid to grids.nc')
            ocpt.oasis.write_grid(
                name=oasis_grid_names['amipfr'],
                lats=amipfr_grid.cell_latitudes(),
                lons=amipfr_grid.cell_longitudes(),
                corners=amipfr_grid.cell_corners()
            )
            self.log_debug('Write AMIP-FR areas to areas.nc')
            ocpt.oasis.write_area(
                name=oasis_grid_names['amipfr'],
                areas=amipfr_grid.cell_areas()
            )
            self.log_debug('Write AMIP-FR mask to masks.nc')
            if self.getarg('amipfr_mask_file', context, default=None):
                self.log_error(
                    'Reading the AMIP-FR mask from file is not implemented yet'
                )
                raise ScriptEngineTaskRunError
            else:
                ocpt.oasis.write_mask(
                    name=oasis_grid_names['amipfr'],
                    masks=np.zeros((amipfr_grid.nx, amipfr_grid.ny))
                )
