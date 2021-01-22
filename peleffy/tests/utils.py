"""
This module contains a variety of helpful tools for tests.
"""


import numpy as np
from simtk import unit

from peleffy.forcefield.forcefield import _BaseForceField


SET_OF_LIGAND_PATHS = ['ligands/benzene.pdb', 'ligands/toluene.pdb',
                       'ligands/oleic_acid.pdb', 'ligands/propionic_acid.pdb',
                       'ligands/trimethylglycine.pdb', 'ligands/ethylene.pdb']


def apply_PELE_dihedral_equation(proper, x):
    """
    Given an x, it applies the PELE's dihedral equation to obtain a y.

    Parameters
    ----------
    proper : a peleffy.topology.Proper
        The proper whose parameters will be applied to equation
    x : float
        Equation's x value
    """
    return proper.constant * (1 + proper.prefactor
                              * np.cos(proper.periodicity * x))


def apply_OFF_dihedral_equation(proper, x):
    """
    Given an x, it applies the Open Force Field's dihedral equation to obtain
    a y.

    Parameters
    ----------
    proper : a peleffy.topology.Proper
        The proper whose parameters will be applied to equation
    x : float
        Equation's x value
    """
    return proper.k * (1 + np.cos(proper.periodicity * x - proper.phase))


def check_CHO_charges(parameters):
    """
    It checks whether C, H and O atoms in a molecule have reasonable
    partial charges assigned.

    Parameters
    ----------
    parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
        The parameter wrapper containing the generated parameters
        along with the partial charges to check

    Raises
    ------
    AssertionError in case that the checking fails
    ValueError if an unexpected element is found in the molecule
    """

    for name, charge in zip(parameters['atom_names'], parameters['charges']):
        charge = charge.value_in_unit(unit.elementary_charge)

        if 'C' in name:
            assert charge < 1.0 and charge > -0.23, \
                'Presumably not a reasonable charge for a C: ' \
                + '{}'.format(charge)
        elif 'H' in name:
            assert charge > 0 and charge < 0.2, \
                'Presumable not a reasonable charge for a H: ' \
                + '{}'.format(charge)
        elif 'O' in name:
            assert charge < -0.5 and charge > -1.0, \
                'Presumable not a reasonable charge for a O: ' \
                + '{}'.format(charge)
        else:
            raise ValueError('Unknown atom name')


def compare_dicts(dict1, dict2):
    """
    Given two dictionaries, it compares them and complains about any
    disagreement.

    Parameters
    ----------
    dict1 : str
        First dictionary to compare
    dict2 : str
        Second dictionary to compare

    Raises
    ------
        AssertionError
            If any difference is found between dictionaries
    """

    assert len(dict1) == len(dict2), 'Number of keys does not match, ' \
        + 'dictionary1: {}, dictionary2: {}'.format(len(dict1), len(dict2))

    for key in dict1.keys():
        assert key in dict2.keys(), 'Key \'{}\' from '.format(key) \
            + 'directory1 not found in dictionary2'

    for key in dict2.keys():
        assert key in dict1.keys(), 'Key \'{}\' from '.format(key) \
            + 'dictionary2 not found in dictionary1'

    for key, value in dict1.items():
        assert value == dict2[key], 'Value for key \'{}\' '.format(key) \
            + 'does not match between dictionaries, ' \
            + 'dictionary1: {}, dictionary2: {}'.format(value, dict2[key])


def merge_dicts(*dict_args):
    """
    Given any number of dictionaries, shallow copy and merge into a new dict,
    precedence goes to key-value pairs in latter dictionaries.

    Parameters
    ----------
    **dict_args : dict
        Dictionary to merge

    Returns
    -------
    merged_dict : str
        Merged dictionary
    """
    merged_dict = {}
    for dictionary in dict_args:
        merged_dict.update(dictionary)
    return merged_dict


def check_parameters(topology, expected_nonbonding=None,
                     expected_bonds=None, expected_angles=None,
                     expected_propers=None, expected_impropers=None):
    """
    It checks the current parameters of the molecule.

    Parameters
    ----------
    topology : a Topology object
        The molecular topology representation to check
    expected_nonbonding : list[list]
        The list of expected nonbonding parameters
    expected_bonds : list[list]
        The list of expected bond parameters
    expected_angles : list[list]
        The list of expected angle parameters
    expected_propers : list[list]
        The list of expected proper parameters
    expected_impropers : list[list]
        The list of expected improper parameters
    """

    from peleffy.template.impact import (WritableAtom,
                                         WritableBond,
                                         WritableAngle,
                                         WritableProper,
                                         WritableImproper)

    if expected_nonbonding is not None:
        assert len(topology.atoms) == len(expected_nonbonding), \
            'Invalid number of nonbonding terms'

        for atom in topology.atoms:
            w_atom = WritableAtom(atom)
            w_parameters = [w_atom.index, w_atom.parent.index, w_atom.core,
                            w_atom.OPLS_type, w_atom.PDB_name,
                            w_atom.unknown, w_atom.sigma, w_atom.epsilon,
                            w_atom.charge, w_atom.born_radius,
                            w_atom.SASA_radius, w_atom.nonpolar_gamma,
                            w_atom.nonpolar_alpha]
            assert w_parameters in expected_nonbonding, \
                'Invalid writable nonbonding parameters ' \
                + '{}'.format(w_parameters)

    if expected_bonds is not None:
        assert len(topology.bonds) == len(expected_bonds), \
            'Invalid number of bond terms'

        for bond in topology.bonds:
            w_bond = WritableBond(bond)
            w_parameters = [attr[1] for attr in list(w_bond)]
            assert w_parameters in expected_bonds, \
                'Invalid writable bond parameters {}'.format(w_parameters)

    if expected_angles is not None:
        assert len(topology.angles) == len(expected_angles), \
            'Invalid number of angle terms'

        for angle in topology.angles:
            w_angle = WritableAngle(angle)
            w_parameters = [attr[1] for attr in list(w_angle)]
            assert w_parameters in expected_angles, \
                'Invalid writable angle parameters {}'.format(w_parameters)

    if expected_propers is not None:
        assert len(topology.propers) == len(expected_propers), \
            'Invalid number of proper terms'

        for proper in topology.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_propers, \
                'Invalid writable proper parameters {}'.format(w_parameters)

    if expected_impropers is not None:
        assert len(topology.impropers) == len(expected_impropers), \
            'Invalid number of improper terms'

        for improper in topology.impropers:
            w_improper = WritableImproper(improper)
            w_parameters = [attr[1] for attr in list(w_improper)]
            assert w_parameters in expected_impropers, \
                'Invalid writable improper parameters {}'.format(w_parameters)


def compare_files(file1, file2):
    """
    It compares two files line by line and gives an AssertionError if there
    is any difference.

    Parameters
    ----------
    file1 : str
        Path to the first file to compare
    file2 : str
        Path to the second file to compare

    Raises
    ------
        AssertionError
            If any difference is found between files
    """
    with open(file1, 'r') as f1:
        lines1 = [line for line in f1.readlines()
                  if not line.startswith("*")]
    with open(file2, 'r') as f2:
        lines2 = [line for line in f2.readlines()
                  if not line.startswith("*")]
    assert len(lines1) == len(lines2), \
        'Number of lines do not match: ' \
        + str(len(lines1)) + ' and ' + str(len(lines2))
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        assert line1 == line2, \
            'Found different lines at line {}:'.format(i) + '\n' + line1 + line2


def compare_files_without_order(file1, file2):
    """
    It compares two files line by line and gives an AssertionError if there
    is any difference. The order of the lines is not taken into account.

    Parameters
    ----------
    file1 : str
        Path to the first file to compare
    file2 : str
        Path to the second file to compare

    Raises
    ------
        AssertionError
            If any difference is found between files
    """
    with open(file1, 'r') as f1:
        lines1 = [line for line in f1.readlines()
                  if not line.startswith("*")]
    with open(file2, 'r') as f2:
        lines2 = [line for line in f2.readlines()
                  if not line.startswith("*")]

    assert len(lines1) == len(lines2), \
        'Number of lines do not match: ' \
        + str(len(lines1)) + ' and ' + str(len(lines2))

    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        assert line1 in lines2, \
            'Line not found in reference template {}:' \
            .format(i) + '\n' + line1


def compare_blocks(block1, block2, line_indices):
    """
    It compares two PDB blocks and gives an AssertionError if there is any
    difference.

    Parameters
    ----------
    block1 : str
        First block to compare
    block2 : str
        Second block to compare
    line_indices : tuple[int, int]
        The line indices to compare

    Raises
    ------
        AssertionError
            If any difference is found between blocks
    """

    idx1, idx2 = line_indices

    lines1 = list()
    lines2 = list()

    for line in block1.split('\n'):
        if line.startswith('REMARK'):
            continue
        if line.startswith('TITLE'):
            continue
        if line.startswith('CONECT'):
            continue
        if line.startswith('CRYST'):
            continue
        if line.startswith('END'):
            continue
        if line.startswith('MODEL'):
            continue
        if line.startswith('ENDMDL'):
            continue
        if line.startswith('TER'):
            continue
        if line == '':
            continue
        lines1.append(line)

    for line in block2.split('\n'):
        if line.startswith('REMARK'):
            continue
        if line.startswith('TITLE'):
            continue
        if line.startswith('CONECT'):
            continue
        if line.startswith('CRYST'):
            continue
        if line.startswith('END'):
            continue
        if line.startswith('MODEL'):
            continue
        if line.startswith('ENDMDL'):
            continue
        if line.startswith('TER'):
            continue
        if line == '':
            continue
        lines2.append(line)

    assert len(lines1) == len(lines2), \
        'Number of lines do not match: ' \
        + str(len(lines1)) + ' and ' + str(len(lines2))

    for line1, line2 in zip(lines1, lines2):
        print(line1[idx1:idx2] + 'X')
        print(line2[idx1:idx2] + 'X')

        assert line1[idx1:idx2] == line2[idx1:idx2], \
            'Found different lines:\n' + line1 + line2


def parameterize_opls2005(opls2005, molecule, ffld_file,
                          charge_method='OPLS2005'):
    """
    This is a workaround to parameterize an OPLS2005 force field
    using a precomputed ffld_file. Thus, we do not require the
    Schrodinger dependency.

    Parameters
    ----------
    opls2005 : an OPLS2005ForceField object
        The OPLS2005 force field to parameterize with the ffld_file
    molecule : a peleffy.topology.Molecule
        The peleffy's Molecule object to parameterize with the ffld file
    ffld_file : str
        The path to the precomputed ffld file from where the parameters
        will be extracted
    charge_method : str
        The charge method to employ to calculate partial charges.
        Default is 'OPLS2005'

    Returns
    -------
    parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
        The parameter wrapper containing the parameters generated
        with the precomputed ffld_file
    """

    from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

    with open(ffld_file) as f:
        ffld_output = f.read()

    parameters = OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                           ffld_output)
    # Assign partial charges using the charge calculator object
    charge_calculator = opls2005._get_charge_calculator(charge_method,
                                                        molecule)
    charge_calculator.assign_partial_charges(parameters)

    return parameters


def parameterize_openffopls2005(openffopls2005, molecule, ffld_file):
    """
    This is a workaround to parameterize an OpenFFOPLS2005 force field
    using a precomputed ffld_file. Thus, we do not require the
    Schrodinger dependency.

    Parameters
    ----------
    openffopls2005 : an OpenFFOPLS2005ForceField object
        The OPLS2005 force field to parameterize with the ffld_file
    molecule : a peleffy.topology.Molecule
        The peleffy's Molecule object to parameterize with the ffld file
    ffld_file : str
        The path to the precomputed ffld file from where the parameters
        will be extracted

    Returns
    -------
    parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
        The parameter wrapper containing the parameters generated
        with the precomputed ffld_file
    """

    parameters = parameterize_opls2005(openffopls2005._oplsff,
                                       molecule,
                                       ffld_file)

    # Initialize mock class
    oplsff = MockOPLS2005ForceField()
    oplsff.set_preloaded_parameters(parameters)

    # Set mock class to the OpenFFOPLS2005ForceField class
    openffopls2005._oplsff = oplsff

    print(parameters['sigmas'])

    return openffopls2005.parameterize(molecule)


class MockBaseForceField(_BaseForceField):
    """
    It is a mock class of _BaseForceField to skip Schrodinger
    dependency.
    """

    _preloaded_parameters = None

    def parameterize(self, molecule):
        """
        It parameterizes the supplied molecule.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """

        # Assign parameters
        parameters = self._preloaded_parameters

        return parameters

    def set_preloaded_parameters(self, parameters):
        """
        It loads parameters to the mock class of BaseForceField.

        Parameters
        ----------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the preloaded parameters
            for this mock class
        """

        self._preloaded_parameters = parameters


class MockOPLS2005ForceField(MockBaseForceField):
    """
    It is a mock class of OPLS2005ForceField to skip Schrodinger
    dependency.
    """

    _preloaded_parameters = None

    def __init__(self):
        """It initializes the OPLS2005 force field class."""
        super().__init__(self._type)

    def _get_parameters(self, molecule):
        """
        Instead of computing the parameters and returning
        """
        return self._preloaded_parameters

    def set_preloaded_parameters(self, parameters):
        """
        It loads parameters to the mock class of OPLS2005ForceField.

        Parameters
        ----------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the preloaded parameters
            for this mock class
        """

        self._preloaded_parameters = parameters
