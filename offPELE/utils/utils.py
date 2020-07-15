from pkg_resources import resource_filename
import os
import contextlib


def get_data_file_path(relative_path):

    fn = resource_filename('offpele', os.path.join(
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


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s:%s' % (filename, lineno, category.__name__, message)


def check_if_path_exists(path):
    from pathlib import Path

    if not isinstance(path, Path):
        path = Path(path)

    if not path.is_dir():
        raise ValueError('Invalid path to {}'.format(path))


def create_path(path):
    import os
    os.makedirs(str(path), exist_ok=True)
