from pkg_resources import resource_filename
import os
import contextlib

from simtk import unit


def get_data_file_path(relative_path):

    fn = resource_filename('offPELE', os.path.join(
        'data', relative_path))

    if not os.path.exists(fn):
        raise ValueError(
            "Sorry! {fn} does not exist. If you just added it, you'll have to "
            + " re-install")

    return fn


@contextlib.contextmanager
def temporary_cd(path):
    old_path = os.getcwd()
    os.chdir(os.path.abspath(path))
    try:
        yield
    finally:
        os.chdir(old_path)


def rmin_halves_to_sigmas(rmin_halves):
    """
    Convert rmin_half values to sigmas according to: http://ambermd.org/Questions/vdwequation.pdf
    """
    FACTOR = 1.122462048309373  # The sixth root of 2

    sigmas = dict()
    for indexes, rmin_half in rmin_halves.items():
        sigma = FACTOR * rmin_half
        sigmas[indexes] = sigma

    return sigmas
