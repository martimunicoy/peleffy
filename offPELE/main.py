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
DEFAULT_OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
DEFAULT_RESOLUTION = int(30)
IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'
DEFAULT_OUTPUT = False
DEFAULT_DATA_LOCAL = False
DEFAULT_SOLVENT = True


# Functions
def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument("pdb_file", metavar="PDB FILE", type=str,
                        help="Path PDB file to parameterize")
    parser.add_argument("-f", "--forcefield", metavar="NAME",
                        type=str, help="OpenForceField's forcefield name. "
                        + "Default is " + str(DEFAULT_OFF_FORCEFIELD),
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


def handle_output_paths(molecule, output, as_datalocal):
    from pathlib import Path
    name = molecule.name
    output_path = Path(output)
    check_if_path_exists(output_path)

    rotlib_path = output_path
    impact_path = output_path
    solvent_path = output_path

    rotlib_name = name.upper() + '.rot.assign'
    impact_name = name.lower() + 'z'
    solvent_name = name.lower() + '_solv.json'

    if as_datalocal:
        rotlib_path = rotlib_path.joinpath(ROTAMER_LIBRARY_PATH)
        impact_path = impact_path.joinpath(IMPACT_TEMPLATE_PATH)
        solvent_path = solvent_path.joinpath(SOLVENT_TEMPLATE_PATH)

        create_path(rotlib_path)
        create_path(impact_path)
        create_path(solvent_path)

    return rotlib_path.joinpath(rotlib_name), \
        impact_path.joinpath(impact_name), \
        solvent_path.joinpath(solvent_name)


def main(pdb_file, forcefield=DEFAULT_OFF_FORCEFIELD, resolution=DEFAULT_RESOLUTION, 
    output=DEFAULT_OUTPUT, as_datalocal=DEFAULT_DATA_LOCAL,
    with_solvent=DEFAULT_SOLVENT):
    print('-' * 60)
    print('Open Force Field parameterizer for PELE v'
          '{}'.format(offPELE.__version__))
    print('-' * 60)
    print(' - PDB to parameterize: {}'.format(pdb_file))
    print(' - Force field: {}'.format(forcefield))
    print(' - Rotamer library resolution: {}'.format(resolution))
    print(' - Output path: {}'.format(output))
    print(' - DataLocal-like output: {}'.format(as_datalocal))
    print('-' * 60)

    # Supress OpenForceField toolkit warnings
    import logging
    logging.getLogger().setLevel(logging.ERROR)

    from offPELE.topology import Molecule
    from offPELE.template import Impact
    from offPELE.solvent import OBC2

    if not output:
        output = os.getcwd()

    molecule = Molecule(pdb_file)
    molecule.parameterize(forcefield)

    rotlib_out, impact_out, solvent_out = handle_output_paths(molecule, output, as_datalocal)

    molecule.build_rotamer_library(resolution=resolution)
    molecule.rotamer_library.to_file(rotlib_out)
    impact = Impact(molecule)
    impact.write(impact_out)

    if with_solvent:
        solvent = OBC2(molecule)
        solvent.to_json_file(solvent_out)

    print(' - All files were generated successfully')
    print('-' * 60)


if __name__ == '__main__':
    args = parse_args()
    main(args.pdb_file, args.forcefield, args.resolution,
        args.output, args.as_datalocal)
