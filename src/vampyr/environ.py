from os import environ
from pathlib import Path
from sysconfig import get_path


def _set_mwfilters_path():
    """
    Sets location of filter files.
    """

    if "MWFILTERS_DIR" not in environ:
        p = Path(get_path("purelib"))
        p = p.parents[2] / "share/MRCPP/mwfilters"
        environ["MWFILTERS_DIR"] = str(p)
