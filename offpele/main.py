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
AVAILABLE_CHARGES_METHODS = ['am1bcc', 'gasteiger', 'OPLS']
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
                        + "compute charges", default=DEFAULT_CHARGES_METHOD,
                        choices=AVAILABLE_CHARGES_METHODS)
    parser.add_argument('--include_terminal_rotamers',
                        dest="include_terminal_rotamers",
                        action='store_true',
                        help="Not exclude terminal rotamers "
                        + "when building the rotamer library")
    parser.add_argument('--use_OPLS_nonbonding_params',
                        dest="use_OPLS_nb_params",
                        action='store_true',
                        help="Use OPLS to set the nonbonding parameters")
    parser.add_argument('--use_OPLS_bonds_and_angles',
                        dest="use_OPLS_bonds_and_angles",
                        action='store_true',
                        help="Use OPLS to set the parameters for bonds "
                        + "and angles")

    parser.set_defaults(as_datalocal=False)
    parser.set_defaults(with_solvent=False)
    parser.set_defaults(include_terminal_rotamers=False)
    parser.set_defaults(use_OPLS_nb_params=False)
    parser.set_defaults(use_OPLS_bonds_and_angles=False)

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
    solvent_name = 'ligandParams.txt'

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
                use_OPLS_nb_params=False,
                use_OPLS_bonds_and_angles=False,
                exclude_terminal_rotamers=True,
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
    use_OPLS_nb_params : bool
        Whether to use Open Force Field or OPLS to obtain the
        nonbonding parameters. Please, note that this option is only
        available if a valid Schrodinger installation is found in the
        current machine. Default is False
    use_OPLS_bonds_and_angles : bool
        Whether to use OPLS to obtain the bond and angle parameters
        or not. Please, note that this option is only
        available if a valid Schrodinger installation is found in the
        current machine. Default is False
    exclude_terminal_rotamers : bool
        Whether to exclude terminal rotamers or not
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
    print('Open Force Field parameterizer for PELE', offpele.__version__)
    print('-' * 60)
    print(' - General:')
    print('   - Input PDB:', pdb_file)
    print('   - Output path:', output)
    print('   - Write solvent parameters:', with_solvent)
    print('   - DataLocal-like output:', as_datalocal)
    print(' - Parameterization:')
    print('   - Force field:', forcefield)
    print('   - Charges method:', charges_method)
    print('   - Use OPLS nonbonding parameters:', use_OPLS_nb_params)
    print('   - Use OPLS bonds and angles:', use_OPLS_bonds_and_angles)
    print(' - Rotamer library:')
    print('   - Resolution:', resolution)
    print('   - Exclude terminal rotamers:', exclude_terminal_rotamers)
    print('-' * 60)

    # Supress OpenForceField toolkit warnings
    import logging
    logging.getLogger().setLevel(logging.ERROR)

    from offpele.topology import Molecule
    from offpele.template import Impact
    from offpele.solvent import OBC2

    if not output:
        output = os.getcwd()

    molecule = Molecule(pdb_file, rotamer_resolution=resolution,
                        exclude_terminal_rotamers=exclude_terminal_rotamers)

    rotlib_out, impact_out, solvent_out = handle_output_paths(molecule,
                                                              output,
                                                              as_datalocal)

    rotamer_library = offpele.topology.RotamerLibrary(molecule)
    rotamer_library.to_file(rotlib_out)

    molecule.parameterize(forcefield, charges_method=charges_method,
                          use_OPLS_nonbonding_params=use_OPLS_nb_params,
                          use_OPLS_bonds_and_angles=use_OPLS_bonds_and_angles)
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

    >>> python main.py molecule.pdb -f openff_unconstrained-1.2.0.offxml
        -r 30 -o output_path/ --with_solvent --as_DataLocal -c gasteiger

    """
    args = parse_args()

    exclude_terminal_rotamers = not args.include_terminal_rotamers

    run_offpele(args.pdb_file, args.forcefield, args.resolution,
                args.charges_method, args.use_OPLS_nb_params,
                args.use_OPLS_bonds_and_angles, exclude_terminal_rotamers,
                args.output, args.with_solvent, args.as_datalocal)


if __name__ == '__main__':
    main()
