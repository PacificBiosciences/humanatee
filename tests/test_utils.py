"""Tests for the humanatee.utils module."""

import pytest
from humanatee import utils

# turn off docstring checking for this file
# ruff: noqa: D103


def test_print_version(capsys):
    utils.print_version('args')
    captured = capsys.readouterr()
    assert 'humanatee v' in captured.out


def test_log_file_stderr(tmp_path, capsys):
    log_file = tmp_path / 'log.txt'
    log_file_stderr = utils.LogFileStderr(log_file)
    log_file_stderr.write('test')
    log_file_stderr.flush()
    assert log_file.read_text() == 'test'
    captured = capsys.readouterr()
    assert 'test' in captured.err


def test_silent_remove():
    assert utils.silent_remove('this_file_doesnt_exist.txt') is None
    assert utils.silent_remove('nonexistent_dir/this_file_doesnt_exist.txt') is None


def test_flatten():
    allele_lengths = [[24, 56], [30, 50], [24, 56], [40], []]
    assert utils.flatten(allele_lengths) == [24, 56, 30, 50, 24, 56, 40]


def test_reverse_counter():
    counter_string = '{"24,56": 2, "30,50": 1, "40": 1}'
    assert utils.reverse_counter(counter_string) == [[24, 56], [24, 56], [30, 50], [40]]


@pytest.mark.parametrize(
    'maybe_int, expected_result',
    [
        ('24', 24),
        ('-56', -56),
        ('0', 0),
        ('', None),
        (' ', None),
        ('a', None),
    ],
)
def test_smart_int(maybe_int, expected_result):
    assert utils.smart_int(maybe_int) == expected_result


@pytest.mark.parametrize(
    'number, test_range, expected_result',
    [
        (24, [0, 100], True),
        (24, [0, 24], True),
        (24, [24, 100], True),
        (24, [25, 100], False),
        (24, [0, 23], False),
        (24, [25, 23], False),
        (24, [24, 24], True),
        (24, [23, 23], False),
    ],
)
def test_in_range(number, test_range, expected_result):
    assert utils.in_range(number, test_range) == expected_result
