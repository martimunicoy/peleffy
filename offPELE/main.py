# -*- coding: utf-8 -*-


# Global imports
import os
import warnings
import argparse as ap

import offPELE
from offPELE.utils import check_if_path_exists, create_path


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Constants
DEFAULT_OFF_FORCEFIELD = 'openff_unconstrained-1.1.1.offxml'
DEFAULT_RESOLUTION = int(30)
IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'


# Functions
def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument("pdb_file", metavar="PDB FILE", type=str,
                        help="Path PDB file to parameterize")
    parser.add_argument("-f", "--forcefield", metavar="OFF FORCEFIELD",
                        type=str, help="OFF forcefield name. Default is "
                        + str(DEFAULT_OFF_FORCEFIELD),
                        default=DEFAULT_OFF_FORCEFIELD)
    parser.add_argument("-r", "--resolution", metavar="INT",
                        help="Rotamer library resolution in degrees. "
                        + "Default is " + str(DEFAULT_RESOLUTION),
                        default=DEFAULT_RESOLUTION)
    parser.add_argument("-o", "--output", metavar="PATH",
                        help="Output path. Default is the current working "
                        + "directory")
    parser.add_argument('--with_solvent', dest='with_solvent',
                        help="Generate solvent parameters for OBC",
                        action='store_true')
    parser.add_argument('--as_DataLocal', dest='as_datalocal',
                        help="Output will be saved following PELE's DataLocal "
                        + "hierarchy", action='store_true')

    parser.set_defaults(as_datalocal=False)
    parser.set_defaults(with_solvent=False)

    args = parser.parse_args()

    if args.output is None:
        args.output = os.getcwd()

    return args


def handle_output_paths(molecule, args):
    from pathlib import Path
    name = molecule.name
    output_path = Path(args.output)
    check_if_path_exists(output_path)

    rotlib_path = output_path
    impact_path = output_path
    solvent_path = output_path

    rotlib_name = name.upper() + '.rot.assign'
    impact_name = name.lower() + 'z'
    solvent_name = name.lower() + '_solv.json'

    if args.as_datalocal:
        rotlib_path = rotlib_path.joinpath(ROTAMER_LIBRARY_PATH)
        impact_path = impact_path.joinpath(IMPACT_TEMPLATE_PATH)
        solvent_path = solvent_path.joinpath(SOLVENT_TEMPLATE_PATH)

        create_path(rotlib_path)
        create_path(impact_path)
        create_path(solvent_path)

    return rotlib_path.joinpath(rotlib_name), \
        impact_path.joinpath(impact_name), \
        solvent_path.joinpath(solvent_name)


def main():
    args = parse_args()
    print('Open Force Field parameterizer for PELE v'
          '{}'.format(offPELE.__version__))
    print('-' * 60)
    print(' - PDB to parameterize: {}'.format(args.pdb_file))
    print(' - Force field: {}'.format(args.forcefield))
    print(' - Rotamer library resolution: {}'.format(args.resolution))
    print(' - Output path: {}'.format(args.output))
    print(' - DataLocal-like output: {}'.format(args.as_datalocal))
    print('-' * 60 + '\n')

    # Supress OpenForceField toolkit warnings
    import logging
    logging.getLogger().setLevel(logging.ERROR)

    from offPELE.topology import Molecule
    from offPELE.template import Impact
    from offPELE.solvent import OBC2

    molecule = Molecule(args.pdb_file)
    molecule.parameterize(args.forcefield)

    rotlib_out, impact_out, solvent_out = handle_output_paths(molecule, args)

    molecule.build_rotamer_library(resolution=args.resolution)
    molecule.rotamer_library.to_file(rotlib_out)
    impact = Impact(molecule)
    impact.write(impact_out)

    if args.with_solvent:
        solvent = OBC2(molecule)
        solvent.to_json_file(solvent_out)

    print(' - All files were generated with successfully')


if __name__ == '__main__':
    main()
