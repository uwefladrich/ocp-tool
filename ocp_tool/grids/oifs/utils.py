import re
import configparser

from collections import namedtuple


def griddes_read(filename):
    """Helper function to provide grid descriptions as returned by 'cdo griddes'
    in a form that can be used to add new grids hardcoded in
    ocp_tool.grids.oifs.
    To use this function, run `cdo griddes OIFS_GRIB_FILE > out_file` and run
    this function on the out_file. The function reads the output and returns a
    dict with one sub-dict for each of the grids. The dictionaries can be
    printed with somthing like
        pp = pprint.PrettyPrinter(width=80, compact=True)
        pp.pprint(d)
    and copy-pasted into new grid description modules."""

    # This is what `cdo -s griddes` is supposed to return:
    griddes_types = dict(
        gridtype  = str,
        gridsize  = int,
        xsize     = int,
        ysize     = int,
        xname     = str,
        xlongname = lambda s: str(s[1:-1]),  # strip quotes
        xunits    = lambda s: str(s[1:-1]),
        yname     = str,
        ylongname = lambda s: str(s[1:-1]),
        yunits    = lambda s: str(s[1:-1]),
        numlpe    = int,
        xvals     = lambda s: [*map(float, s.split())],  # list of floats
        yvals     = lambda s: [*map(float, s.split())],
        reducedpoints = lambda s: [*map(int, s.split())],  # list of ints
        rowlon    = lambda s: [*map(int, s.split())],  # list of ints
    )

    cfg = configparser.ConfigParser()
    with open(filename) as f:
        cfg.read_string(
            re.sub(
                r'#\s*\n#\s+(.*)\n#\s*\n',  # Matches the cdo griddes headers
                r'[\1]\n',             # and replaces by configparser headers
                f.read(),
                re.MULTILINE
            )
        )

    griddes = dict()
    for sec in cfg.sections():
        griddes[sec] = dict()
        for key, val in cfg[sec].items():
            try:
                griddes[sec][key] = griddes_types[key](val)
            except KeyError:
                griddes[sec][key] = val

    return griddes


def namedtuple_from_dict(name, dict_):
    return namedtuple(name, dict_.keys())(**dict_)
