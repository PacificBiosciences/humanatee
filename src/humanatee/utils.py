"""Miscellaneous utilities."""

import contextlib
import errno
import gzip
import json
import logging
import os
import sys
import warnings
from collections import Counter
from importlib.metadata import version


def print_version(args):
    """Print the version."""
    sys.stdout.write(f"humanatee v{version('humanatee')}\n")


class LogFileStderr:
    """Write to stderr and a file.

    Useful in conjunction with :meth:`setup_logging`

    Example: setup_logging(stream=LogFileStderr('log.txt'))
    """

    def __init__(self, filename: str):
        """Keep these props."""
        self.name = filename
        self.file_handler = open(filename, 'w')

    def write(self, *args: str) -> None:
        """Write to handlers."""
        sys.stderr.write(*args)
        self.file_handler.write(*args)

    def flush(self) -> None:
        """Flush handlers."""
        sys.stderr.flush()
        self.file_handler.flush()


def setup_logging(
    debug=False,
    stream=sys.stderr,
    log_format='%(asctime)s [%(levelname)s] %(message)s',
    show_version=False,
) -> None:
    """Create default logger."""
    logLevel = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    if show_version:
        logging.info(f"humanatee v{version('humanatee')}")
        logging.info(f"Command: {' '.join(sys.argv)}")

    def sendWarningsToLog(message, category, filename, lineno, *args, **kwargs) -> None:
        """Put warnings into logger."""
        logging.warning(f'{filename}:{lineno}: {category.__name__}:{message}')

    warnings.showwarning = sendWarningsToLog


def writer(fn):
    """Write to file or stdout."""

    @contextlib.contextmanager
    def stdout():
        yield sys.stdout

    if fn:
        if fn.lower().endswith('.gz'):
            return gzip.open(fn, 'wt')
        else:
            return open(fn, 'w')
    else:
        return stdout()


def silent_remove(filename):
    """Remove a file and ignore if it doesn't exist."""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


def flatten(nested_list):
    """Flatten a nested list."""
    return [x for xs in nested_list for x in xs]


def reverse_counter(counter_string):
    """Revert a counter string to a list of elements."""
    return [
        list(map(smart_int, genotype.split(','))) for genotype in list(Counter(json.loads(counter_string)).elements())
    ]


def smart_int(x):
    """Convert to int if possible."""
    try:
        return int(x)
    except ValueError:
        return None


def in_range(x, range):
    """Check if x is in range."""
    return range[0] <= x <= range[1]
