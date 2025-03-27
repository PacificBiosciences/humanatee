"""Main entrypoint for humanatee tools."""

import argparse
import sys
from importlib.metadata import version

from humanatee import print_version, trgt_main, trgtpop_main


def main(prog='humanatee', level=1):
    """Run from command line."""
    CMDS = {
        'trgt': ('Annotate TRGT results', trgt_main),
        'trgtpop': ('Aggregate population VCFs for annotation', trgtpop_main),
        'version': ('Print the version and exit', print_version),
    }

    USAGE = '\n'.join(
        [
            f"humanatee v{version('humanatee')} HiFi Human Annotation Tools",
            '\n',
            'commands:',
            '\n'.join([f'  {k:11} {t[0]}' for k, t in CMDS.items()]),
        ]
    )

    parser = argparse.ArgumentParser(
        prog=prog,
        description=USAGE,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        'cmd',
        metavar='CMD',
        choices=CMDS.keys(),
        type=str,
        default=None,
        help='Command to execute',
    )
    parser.add_argument(
        'options',
        metavar='OPTIONS',
        nargs=argparse.REMAINDER,
        help='Options to pass to the command',
    )

    if len(sys.argv) == level:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args(sys.argv[level:])

    CMDS[args.cmd][1](args.options)


if __name__ == '__main__':
    main()
