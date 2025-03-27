"""Miscellaneous DButilities."""

import logging
import sqlite3

import pandas as pd


def smart_execute(cursor, query, values=(), many=False, always_print=False):
    """Print SQL query if execution fails."""
    if always_print:
        logging.info(f'{query}\n{values}')
    if many:
        try:
            cursor.executemany(query, values)
        except Exception as e:
            logging.error(f'Error executing query:\n{query}')
            raise e
    else:
        try:
            cursor.execute(query, values)
        except Exception as e:
            logging.error(f'Error executing query:\n{query}')
            raise e
    return cursor


def validate_db_schema(conn1, cursor1, conn2, db2):
    """Validate that db2 has the same tables and table info as db in conn1."""
    cursor1.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor1.fetchall()
    for table in tables:
        expected = pd.read_sql_query(f'pragma table_info({table[0]})', conn1)
        actual = pd.read_sql_query(f'pragma table_info({table[0]})', conn2)
        if not expected.equals(actual):
            logging.error(f"Table '{table[0]}' in database '{db2}' does not match expected table definition.\n")
            logging.debug(f'Expected:\n{expected}\n')
            logging.debug(f'Actual:\n{actual}\n')
            raise ValueError(f"Table '{table[0]}' in database '{db2}' does not match expected table definition.")


def add_column(cursor, table, column, type):
    """Add column to table if it does not already exist."""
    try:
        cursor.execute(f'ALTER TABLE {table} ADD COLUMN {column} {type};')
    except sqlite3.OperationalError as e:
        if 'duplicate column name' in str(e):
            pass
        else:
            raise e


def add_table_row(cursor, table, fields_dict):
    """Add a row to a SQL table with fields defined in dict."""
    if fields_dict == {}:
        return
    cursor.execute(
        f"""
        INSERT INTO
        {table} ({', '.join(fields_dict.keys())})
        VALUES (? {', ?' * (len(fields_dict) - 1)})
        """,
        (*fields_dict.values(),),
    )
