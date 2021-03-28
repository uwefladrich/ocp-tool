from .gaussian import GaussianGrid, ReducedGaussianGrid
from .orca import ORCA
from .ifs_griddes import (
    TCO95, O96,
    TCO159, O160,
    TL159, N80,
    TL255, N128,
)


def factory(type_, *args, **kwargs):

    full_gaussian_grids = {
        'F128': N128,
    }

    if type_ in (
        TCO95, O96,
        TCO159, O160,
        TL159, N80,
        TL255, N128,
    ):
        return ReducedGaussianGrid(
            nlats=type_.ysize,
            lats=type_.yvals,
            nlons=type_.reducedpoints
        )
    elif type_ in ('ORCA', 'orca'):
        return ORCA(*args, *kwargs)

    elif type_ in full_gaussian_grids:
        return GaussianGrid(
            lats=full_gaussian_grids[type_].yvals
        )

    raise NotImplementedError(f'Unknown grid type: {type_}')
