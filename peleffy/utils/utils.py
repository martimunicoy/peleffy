"""
This module contains miscellaneous set of handy classes and functions.
"""


__all__ = ["get_data_file_path",
           "temporary_cd",
           "warning_on_one_line",
           "check_if_path_exists",
           "create_path",
           "Logger"
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
    output_path = resource_filename('peleffy', os.path.join(
        'data', relative_path))

    if not os.path.exists(output_path):
        raise ValueError(
            "Sorry! {} does not exist. ".format(output_path)
            + "If you just added it, you'll have to re-install")

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


class Logger(object):
    """
    It contains all the required methods to handle logging messages.
    """
    import logging
    DEFAULT_LEVEL = logging.INFO

    def __init__(self):
        """It initializes a Logger object"""
        import logging

        # Get peleffy logger and set level only the first time
        if 'peleffy_log' not in logging.root.manager.loggerDict:
            self._logger = logging.getLogger('peleffy_log')
            self._logger.setLevel(self.DEFAULT_LEVEL)
        else:
            self._logger = logging.getLogger('peleffy_log')

        # If no handler is found add stream handler
        if not len(self._logger.handlers):
            ch = logging.StreamHandler()
            ch.setLevel(self.DEFAULT_LEVEL)
            formatter = logging.Formatter('%(message)s')
            ch.setFormatter(formatter)
            self._logger.addHandler(ch)

    def set_level(self, level):
        """
        It sets the logging level.

        Parameters
        ----------
        level : str
            The logging level to set. One of [DEBUG, INFO, WARNING, ERROR,
            CRITICAL]
        """
        import logging

        if level.upper() == 'DEBUG':
            logging_level = logging.DEBUG
        elif level.upper() == 'INFO':
            logging_level = logging.INFO
        elif level.upper() == 'WARNING':
            logging_level = logging.WARNING
        elif level.upper() == 'ERROR':
            logging_level = logging.ERROR
        elif level.upper() == 'CRITICAL':
            logging_level = logging.CRITICAL
        else:
            raise ValueError('Invalid level type')

        self._logger.setLevel(logging_level)
        for handler in self._logger.handlers:
            handler.setLevel(logging_level)

    def debug(self, *messages):
        """
        It pulls a debug message.

        Parameters
        ----------
        messages : list[str]
            The list of messages to print
        """
        if len(messages) > 1:
            self._logger.debug(' '.join(map(str, messages)))
        else:
            self._logger.debug(messages[0])

    def info(self, *messages):
        """
        It pulls an info message.

        Parameters
        ----------
        messages : list[str]
            The list of messages to print
        """
        if len(messages) > 1:
            self._logger.info(' '.join(map(str, messages)))
        else:
            self._logger.info(messages[0])

    def warning(self, *messages):
        """
        It pulls a warning message.

        Parameters
        ----------
        messages : list[str]
            The list of messages to print
        """
        if len(messages) > 1:
            self._logger.warning(' '.join(map(str, messages)))
        else:
            self._logger.warning(messages[0])

    def error(self, *messages):
        """
        It pulls a error message.

        Parameters
        ----------
        messages : list[str]
            The list of messages to print
        """
        if len(messages) > 1:
            self._logger.error(' '.join(map(str, messages)))
        else:
            self._logger.error(messages[0])

    def critical(self, *messages):
        """
        It pulls a critical message.

        Parameters
        ----------
        messages : list[str]
            The list of messages to print
        """
        if len(messages) > 1:
            self._logger.critical(' '.join(map(str, messages)))
        else:
            self._logger.critical(messages[0])
