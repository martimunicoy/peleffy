from pkg_resources import resource_filename
import os
import contextlib


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
