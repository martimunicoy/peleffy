# -*- coding: utf-8 -*-


# Global imports
import argparse as ap


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument("pdb_file", metavar="FILE", type=str,
                        help="Path PDB file to parameterize")
    args = parser.parse_args()

    return args.pdb_file


def main():
    pass


if __name__ == '__main__':
    main()
