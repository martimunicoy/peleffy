# -*- coding: utf-8 -*-
"""
This module is designed to run offpele through the command-line.
"""

__author__ = "Marti Municoy"
__license__ = "GPL"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


import os
import argparse as ap

import offpele
from offpele.utils import check_if_path_exists, create_path


DEFAULT_OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
DEFAULT_RESOLUTION = int(30)
DEFAULT_CHARGES_METHOD = 'am1bcc'
IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'
DEFAULT_TERMINAL_ROT_TO_IGNORE = 1


def parse_args():
    """
    It parses the command-line arguments.

    Returns
    -------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """
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
    parser.add_argument('-c', '--charges_method', metavar="NAME",
                        type=str, help="The name of the method to use to "
                        + "compute charges", default=DEFAULT_CHARGES_METHOD)
    parser.add_argument('-t', '--terminal_rotamers_to_ignore', metavar="INT",
                        type=str, help="The number of terminal rotamers " +
                        " to ignore when building the rotamer library",
                        default=DEFAULT_TERMINAL_ROT_TO_IGNORE)

    parser.set_defaults(as_datalocal=False)
    parser.set_defaults(with_solvent=False)

    args = parser.parse_args()

    return args


def handle_output_paths(molecule, output, as_datalocal):
    """
    It handles the output paths where offpele's output files will be saved.

    Parameters
    ----------
    molecule : offpele.topology.Molecule
        A Molecule object
    output : str
        The output path supplied by the user
    as_datalocal : bool
        Whether to save output files following PELE's DataLocal hierarchy or
        not

    Returns
    -------
    rotlib_path : pathlib.Path
        The output path for the rotamer library
    impact_path : pathlib.Path
        The output path for the Impact template
    solvent_path : pathlib.Path
        The output path for the solvent template
    """
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


def run_offpele(pdb_file, forcefield=DEFAULT_OFF_FORCEFIELD,
                resolution=DEFAULT_RESOLUTION,
                charges_method=DEFAULT_CHARGES_METHOD,
                terminal_rotamers_to_ignore=DEFAULT_TERMINAL_ROT_TO_IGNORE,
                output=None, with_solvent=False, as_datalocal=False):
    """
    It runs offpele.

    Parameters
    ----------
    pdb_file : str
        The path to the pdb_file to parameterize with offpele
    forcefield : str
        The name of an OpenForceField's forcefield
    resolution : float
        The resolution in degrees for the rotamer library. Default is 30
    charges_method : str
        The name of the method to use to compute partial charges. Default
        is 'am1bcc'
    terminal_rotamers_to_ignore : int
        The number of terminal rotamers to ignore when building the
        rotamer library. Default is 1
    output : str
        Path where output files will be saved
    with_solvent : bool
        Whether to generate and save the solvent parameters for the input
        molecule or not
    as_datalocal : bool
        Whether to save output files following PELE's DataLocal hierarchy or
        not
    """
    print('-' * 60)
    print('Open Force Field parameterizer for PELE '
          '{}'.format(offpele.__version__))
    print('-' * 60)
    print(' - PDB to parameterize: {}'.format(pdb_file))
    print(' - Force field: {}'.format(forcefield))
    print(' - Rotamer library resolution: {}'.format(resolution))
    print(' - Charges method: {}'.format(charges_method))
    print(' - Terminal rotamers to ignore: {}'.format(
        terminal_rotamers_to_ignore))
    print(' - Output path: {}'.format(output))
    print(' - Write solvent parameters: {}'.format(with_solvent))
    print(' - DataLocal-like output: {}'.format(as_datalocal))
    print('-' * 60)

    # Supress OpenForceField toolkit warnings
    import logging
    logging.getLogger().setLevel(logging.ERROR)

    from offpele.topology import Molecule
    from offpele.template import Impact
    from offpele.solvent import OBC2

    if not output:
        output = os.getcwd()

    molecule = Molecule(pdb_file)
    molecule.parameterize(forcefield, charges_method=charges_method)

    rotlib_out, impact_out, solvent_out = handle_output_paths(molecule, output, as_datalocal)

    molecule.build_rotamer_library(
        resolution=resolution,
        n_rot_bonds_to_ignore=terminal_rotamers_to_ignore)
    molecule.rotamer_library.to_file(rotlib_out)
    impact = Impact(molecule)
    impact.write(impact_out)

    if with_solvent:
        solvent = OBC2(molecule)
        solvent.to_json_file(solvent_out)

    print(' - All files were generated successfully')
    print('-' * 60)


def main():
    """
    It reads the command-line arguments and runs offpele.

    Examples
    --------

    From the command-line:

    >>> python main.py molecule.pdb -f openff_unconstrained-1.1.1.offxml -r 30
        -o output_path/ --with_solvent --as_DataLocal -c gasteiger

    """
    args = parse_args()
    run_offpele(args.pdb_file, args.forcefield, args.resolution,
                args.charges_method, args.terminal_rotamers_to_ignore,
                args.output, args.with_solvent,
                args.as_datalocal)


if __name__ == '__main__':
    main()
