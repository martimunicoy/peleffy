"""
This module contains miscellaneous set of handy classes and functions.
"""


__all__ = ["get_data_file_path",
           "temporary_cd",
           "warning_on_one_line",
           "check_if_path_exists",
           "create_path"
           ]


from pkg_resources import resource_filename
import os
import contextlib


def get_data_file_path(relative_path):
    """
    It returns the path in the package's data location.

    Parameters
    ----------
    relative_path : str
        The relative path to the file that is required

    Returns
    -------
    output_path : str
        The path in the package's data location, if found
    """
    output_path = resource_filename('offpele', os.path.join(
        'data', relative_path))

    if not os.path.exists(output_path):
        raise ValueError(
            "Sorry! {output_path} does not exist. If you just added it, "
            + "you'll have to re-install")

    return output_path


@contextlib.contextmanager
def temporary_cd(path):
    """
    A context that applies a temporary move to certain path.

    Parameters
    ----------
    path : str
        The path to move to
    """
    old_path = os.getcwd()
    os.chdir(os.path.abspath(path))
    try:
        yield
    finally:
        os.chdir(old_path)


def warning_on_one_line(message, category, filename, lineno, file=None,
                        line=None):
    """A customized warning output in a single line."""
    return ' %s:%s: %s:%s' % (filename, lineno, category.__name__, message)


def check_if_path_exists(path):
    """
    It checks if the supplied path exists and raises a ValueError otherwise.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to check
    """
    from pathlib import Path

    if not isinstance(path, Path):
        path = Path(path)

    if not path.is_dir():
        raise ValueError('Invalid path to {}'.format(path))


def create_path(path):
    """
    It creates a path.

    Parameters
    ----------
    path : str
        The path that will be created
    """
    os.makedirs(str(path), exist_ok=True)
