"""Miscellaneous DButilities."""

import logging
import os
import sqlite3
from importlib import resources

from humanatee import utils


def smart_execute(cursor, query, values=(), many=False, print_query='error'):
    """Print SQL query if execution fails."""
    if print_query is True:
        logging.info(f'{query}\n{values}')
    if many:
        try:
            cursor.executemany(query, values)
        except Exception as e:
            logging.error(f'Error executing query:\n{query}') if print_query is not False else None
            raise e
    else:
        try:
            cursor.execute(query, values)
        except Exception as e:
            logging.error(f'Error executing query:\n{query}') if print_query is not False else None
            raise e
    return cursor


def validate_db_tables(cursor, table_dict, raise_error=True):
    """Validate that db has expected tables and columns."""
    tables_valid = True
    columns_valid = True
    attached = ''
    for table, required_columns in table_dict.items():
        if '.' in table:
            attached, table = table.split('.')
            attached = f'{attached}.'
        db_tables = cursor.execute(f"SELECT name FROM {attached}sqlite_master WHERE type='table';").fetchall()
        if table not in [table[0] for table in db_tables]:
            tables_valid = False
            if raise_error:
                raise ValueError(f'Database does not have required table ({table}).')
        columns = cursor.execute(f'PRAGMA {attached}table_info("{table}");').fetchall()
        columns = [column[1] for column in columns]
        if not all([col in columns for col in required_columns]):
            columns_valid = False
            if raise_error:
                raise ValueError(f'Table ({table}) does not have required columns.')
    return tables_valid and columns_valid


def add_columns(cursor, table, column_type_dict, raise_error=False):
    """Add column to table if it does not already exist."""
    for column, sql_type in column_type_dict.items():
        try:
            cursor.execute(f'ALTER TABLE {table} ADD COLUMN {column} {sql_type};')
        except sqlite3.OperationalError as e:
            if 'duplicate column name' in str(e) and not raise_error:
                pass
            else:
                raise e


def add_table_row(cursor, table, fields_dict, update_conflict='', replace=False, print_query='error'):
    """Add a row to a SQL table with fields defined in dict."""
    if fields_dict == {}:
        return
    if update_conflict and replace:
        raise ValueError('Cannot use both update and replace options.')
    if update_conflict:
        update_conflict = f"""
            ON CONFLICT({update_conflict}) DO UPDATE SET
            {', '.join([f'{k}=excluded.{k}' for k in fields_dict.keys() if k != {update_conflict}])}
            """
    command = 'INSERT OR REPLACE' if replace else 'INSERT'
    smart_execute(
        cursor,
        f"""
        {command} INTO
        {table} ({', '.join(fields_dict.keys())})
        VALUES (? {', ?' * (len(fields_dict) - 1)})
        {update_conflict}
        """,
        (*fields_dict.values(),),
        print_query=print_query,
    )


def initialize_db(prefix, schemas=[], delete=True, no_db=False):
    """Initialize a SQLite database with the given schema."""
    # initialize output database, removing if it already exists
    logging.debug(f'Initializing database with prefix: {prefix}')
    db_path = ''  # see https://stackoverflow.com/a/32833770 instead of ':memory:'
    if not no_db:
        db_path = f'{prefix}.db'
        # create parent directories if necessary
        if os.path.dirname(db_path) != '':
            os.makedirs(os.path.dirname(db_path), exist_ok=True)
    if delete:
        utils.silent_remove(db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    for schema in schemas:
        vcf_schema = resources.files('humanatee').joinpath('schemas', f'{schema}.sql').read_text()
        cursor.executescript(vcf_schema)
    return conn, cursor
