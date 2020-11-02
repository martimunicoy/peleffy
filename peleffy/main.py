# -*- coding: utf-8 -*-
"""
This module is designed to run peleffy through the command-line.
"""

__author__ = "Marti Municoy"
__license__ = "GPL"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


import os
import argparse as ap

import peleffy
from peleffy.utils import Logger, OutputPathHandler


DEFAULT_OFF_FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'
DEFAULT_RESOLUTION = int(30)
DEFAULT_CHARGE_METHOD = 'am1bcc'
AVAILABLE_CHARGE_METHODS = ['am1bcc', 'gasteiger', 'OPLS']
IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'
DEFAULT_TERMINAL_ROT_TO_IGNORE = 1


def parse_args(args):
    """
    It parses the command-line arguments.

    Parameters
    ----------
    args : list[str]
        List of command-line arguments to parse

    Returns
    -------
    parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """
    parser = ap.ArgumentParser()
    parser.add_argument("pdb_file", metavar="PDB FILE", type=str,
                        help="Path PDB file to parameterize")
    parser.add_argument("-f", "--forcefield", metavar="NAME",
                        type=str, help="OpenForceField's forcefield name. "
                        + "Default is " + str(DEFAULT_OFF_FORCEFIELD),
                        default=DEFAULT_OFF_FORCEFIELD)
    parser.add_argument("-r", "--resolution", metavar="INT", type=int,
                        help="Rotamer library resolution in degrees. "
                        + "Default is " + str(DEFAULT_RESOLUTION),
                        default=DEFAULT_RESOLUTION)
    parser.add_argument("-o", "--output", metavar="PATH",
                        help="Output path. Default is the current working "
                        + "directory")
    parser.add_argument('--with_solvent', dest='with_solvent',
                        help="Generate solvent parameters for OBC",
                        action='store_true')
    parser.add_argument('--as_datalocal', dest='as_datalocal',
                        help="Output will be saved following PELE's DataLocal "
                        + "hierarchy", action='store_true')
    parser.add_argument('-c', '--charge_method', metavar="NAME",
                        type=str, help="The name of the method to use to "
                        + "compute charges", default=DEFAULT_CHARGE_METHOD,
                        choices=AVAILABLE_CHARGE_METHODS)
    parser.add_argument('--include_terminal_rotamers',
                        dest="include_terminal_rotamers",
                        action='store_true',
                        help="Not exclude terminal rotamers "
                        + "when building the rotamer library")
    parser.add_argument('-s', '--silent',
                        dest="silent",
                        action='store_true',
                        help="Activate silent mode")
    parser.add_argument('-d', '--debug',
                        dest="debug",
                        action='store_true',
                        help="Activate debug mode")

    parser.set_defaults(as_datalocal=False)
    parser.set_defaults(with_solvent=False)
    parser.set_defaults(include_terminal_rotamers=False)
    parser.set_defaults(silent=False)
    parser.set_defaults(debug=False)

    parsed_args = parser.parse_args(args)

    return parsed_args


def run_peleffy(pdb_file, forcefield=DEFAULT_OFF_FORCEFIELD,
                resolution=DEFAULT_RESOLUTION,
                charge_method=DEFAULT_CHARGE_METHOD,
                exclude_terminal_rotamers=True,
                output=None, with_solvent=False, as_datalocal=False):
    """
    It runs peleffy.

    Parameters
    ----------
    pdb_file : str
        The path to the pdb_file to parameterize with peleffy
    forcefield : str
        The name of an OpenForceField's forcefield
    resolution : float
        The resolution in degrees for the rotamer library. Default is 30
    charge_method : str
        The name of the method to use to compute partial charges. Default
        is 'am1bcc'
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
    log = Logger()
    log.info('-' * 60)
    log.info('Open Force Field parameterizer for PELE', peleffy.__version__)
    log.info('-' * 60)
    log.info(' - General:')
    log.info('   - Input PDB:', pdb_file)
    log.info('   - Output path:', output)
    log.info('   - Write solvent parameters:', with_solvent)
    log.info('   - DataLocal-like output:', as_datalocal)
    log.info(' - Parameterization:')
    log.info('   - Force field:', forcefield)
    log.info('   - Charge method:', charge_method)
    log.info(' - Rotamer library:')
    log.info('   - Resolution:', resolution)
    log.info('   - Exclude terminal rotamers:', exclude_terminal_rotamers)
    log.info('-' * 60)

    from peleffy.topology import Molecule
    from peleffy.template import Impact
    from peleffy.solvent import OBC2

    if not output:
        output = os.getcwd()

    molecule = Molecule(pdb_file, rotamer_resolution=resolution,
                        exclude_terminal_rotamers=exclude_terminal_rotamers)

    output_handler = OutputPathHandler(molecule, output_path=output,
                                       as_datalocal=as_datalocal)

    rotamer_library = peleffy.topology.RotamerLibrary(molecule)
    rotamer_library.to_file(output_handler.get_rotamer_library_path())

    molecule.parameterize(forcefield, charge_method=charge_method)
    impact = Impact(molecule)
    impact.write(output_handler.get_impact_template_path())

    if with_solvent:
        solvent = OBC2(molecule)
        solvent.to_json_file(output_handler.get_solvent_template_path())

    log.info(' - All files were generated successfully')
    log.info('-' * 60)


def main(args):
    """
    It reads the command-line arguments and runs peleffy.

    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user

    Examples
    --------

    From the command-line:

    >>> python main.py molecule.pdb -f openff_unconstrained-1.2.0.offxml
        -r 30 -o output_path/ --with_solvent --as_datalocal -c gasteiger

    """

    exclude_terminal_rotamers = not args.include_terminal_rotamers

    # Supress OpenForceField toolkit warnings
    import logging
    logging.getLogger().setLevel(logging.ERROR)

    # Set peleffy logger to the corresponding level
    logger = Logger()
    if args.silent:
        logger.set_level('CRITICAL')
    elif args.debug:
        logger.set_level('DEBUG')
    else:
        logger.set_level('INFO')

    run_peleffy(args.pdb_file, args.forcefield, args.resolution,
                args.charge_method, exclude_terminal_rotamers,
                args.output, args.with_solvent, args.as_datalocal)


if __name__ == '__main__':
    import sys
    args = parse_args(sys.argv[1:])
    main(args)
