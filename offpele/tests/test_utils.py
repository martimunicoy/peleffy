"""
This module contains the tests to check some handy classes and functions
of offpele.
"""


import pytest

import io
from contextlib import redirect_stdout
from offpele.utils import Logger


class TestLogger(object):
    def test_logger_levels(self):
        """
        It checks the correct behaviour of the different log levels.
        """
        def push_messages(log):
            """Pull some messages at different levels."""
            log.debug('Debug message')
            log.info('Info message')
            log.warning('Warn message')
            log.error('Error message')
            log.critical('Critical message')

        import logging

        # Initiate logger
        log = Logger()

        # Try the default level (INFO)
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Info message\nWarn message\n' \
                + 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try DEBUG level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try DEBUG level
            log.set_level('DEBUG')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Debug message\nInfo message\n'\
                + 'Warn message\nError message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try INFO level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try INFO level
            log.set_level('INFO')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Info message\nWarn message\n' \
                + 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try WARNING level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try WARNING level
            log.set_level('WARNING')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Warn message\nError message\n' \
                + 'Critical message\n', \
                'Unexpected logger message at standard output'

        # Try ERROR level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try ERROR level
            log.set_level('ERROR')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try CRITICAL level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try CRITICAL level
            log.set_level('CRITICAL')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Critical message\n', \
                'Unexpected logger message at standard output'
