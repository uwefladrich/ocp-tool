from .gaussian import ReducedGaussianGrid
from .orca import ORCA
from .ifs_griddes import (
    TQ21, F16,
    TCO95, O96,
    TCO159, O160,
    TL159, N80,
    TL255, N128,
)


def factory(type_, *args, **kwargs):

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

    raise NotImplementedError(f'Unknown grid type: {type_}')
