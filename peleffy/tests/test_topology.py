"""
This module contains the tests to check all available force fields in
peleffy.
"""


class TestTopology(object):
    """
    It wraps all tests that check the OpenForceField class.
    """

    def test_empty_parameters(self):
        """
        It tests the initialization of a Topology object with an empty
        parameters wrapper.
        """

        from peleffy.topology import Molecule
        from peleffy.forcefield.parameters import BaseParameterWrapper
        from peleffy.topology import Topology

        molecule = Molecule()
        parameters = BaseParameterWrapper()
        Topology(molecule, parameters)

    def test_add_topological_elements(self):
        """
        It tests the addition of topological elements to an empty
        topology.
        """

        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField
        from peleffy.forcefield.parameters import BaseParameterWrapper
        from peleffy.topology import Topology
        from peleffy.utils import get_data_file_path

        # Define molecule1 and its topology
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        molecule1 = Molecule(pdb_path)
        openff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        parameters1 = openff.parameterize(molecule1)
        topology1 = Topology(molecule1, parameters1)

        # Define empty topology2
        molecule2 = Molecule()
        parameters2 = BaseParameterWrapper()
        topology2 = Topology(molecule2, parameters2)

        # Add parameters to topology2
        for atom in topology1.atoms:
            topology2.add_atom(atom)
        for bond in topology1.bonds:
            topology2.add_bond(bond)
        for angle in topology1.angles:
            topology2.add_angle(angle)
        for proper in topology1.propers:
            topology2.add_proper(proper)
        for improper in topology1.impropers:
            topology2.add_improper(improper)

        # Verify content of both topologies
        assert topology1.atoms == topology2.atoms, \
            'The atoms of boths topologies should match'
        assert topology1.bonds == topology2.bonds, \
            'The bonds of boths topologies should match'
        assert topology1.angles == topology2.angles, \
            'The angles of boths topologies should match'
        assert topology1.propers == topology2.propers, \
            'The propers of boths topologies should match'
        assert topology1.impropers == topology2.impropers, \
            'The impropers of boths topologies should match'

    def test_openff_parameterizer(self):
        """
        It checks the behaviour of the Topology with the OpenFF
        parameters.
        """

        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField
        from peleffy.topology import Topology
        from peleffy.utils import get_data_file_path
        from .utils import check_parameters

        FORCE_FIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

        # Load molecule
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        openff = OpenForceField(FORCE_FIELD_NAME)

        # Obtain force field parameters
        parameters = openff.parameterize(molecule)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_nonbonding = [
            [1, 0, 'M', 'OFFT', '_C1_', 0, 3.3996695084235347, 0.1094,
             -0.1088, 0, 1.6998347542117673, 0, 0],
            [2, 1, 'M', 'OFFT', '_H1_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [3, 1, 'M', 'OFFT', '_H2_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [4, 1, 'M', 'OFFT', '_H3_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [5, 1, 'M', 'OFFT', '_H4_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0]]

        expected_bonds = [
            [1, 2, 376.8940758588, 1.094223427522],
            [1, 3, 376.8940758588, 1.094223427522],
            [1, 4, 376.8940758588, 1.094223427522],
            [1, 5, 376.8940758588, 1.094223427522]]

        expected_angles = [
            [2, 1, 3, 33.78875634641, 110.2468561538],
            [2, 1, 4, 33.78875634641, 110.2468561538],
            [2, 1, 5, 33.78875634641, 110.2468561538],
            [3, 1, 4, 33.78875634641, 110.2468561538],
            [3, 1, 5, 33.78875634641, 110.2468561538],
            [4, 1, 5, 33.78875634641, 110.2468561538]]

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_nonbonding,
                         expected_bonds=expected_bonds,
                         expected_angles=expected_angles)

        # Load molecule
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        openff = OpenForceField(FORCE_FIELD_NAME)

        # Obtain force field parameters
        parameters = openff.parameterize(molecule)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_propers = [
            [3, 1, 2, 5, 5.376019778605, -1, 2, 0.0],
            [3, 1, 2, 6, 5.376019778605, -1, 2, 0.0],
            [4, 1, 2, 5, 5.376019778605, -1, 2, 0.0],
            [4, 1, 2, 6, 5.376019778605, -1, 2, 0.0]]

        expected_impropers = [
            [1, 2, 5, 6, 1.1, -1, 2],
            [2, 1, 3, 4, 1.1, -1, 2]]

        # Check it up
        check_parameters(topology,
                         expected_propers=expected_propers,
                         expected_impropers=expected_impropers)

    def test_opls2005_parameterizer(self):
        """
        It checks the behaviour of the Topology with the OPLS2005
        parameters.
        """

        from peleffy.topology import Molecule
        from peleffy.topology import Topology
        from peleffy.forcefield import OPLS2005ForceField
        from peleffy.utils import get_data_file_path
        from .utils import check_parameters, parameterize_opls2005

        # Load molecule
        opls2005 = OPLS2005ForceField()
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        ffld_file = get_data_file_path('tests/MET_ffld_output.txt')

        # Workaround to parameterize with OPLS2005 without the Schrodinger
        # dependency
        parameters = parameterize_opls2005(opls2005, molecule, ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_nonbonding = [
            [1, 0, 'M', 'CT', '_C1_', 0, 3.5, 0.066, -0.24, 1.975, 1.75,
             0.005, -0.74168571],
            [2, 1, 'M', 'HC', '_H1_', 0, 2.5, 0.03, 0.06, 1.425, 1.25,
             0.00859824, 0.268726247],
            [3, 1, 'M', 'HC', '_H2_', 0, 2.5, 0.03, 0.06, 1.425, 1.25,
             0.00859824, 0.268726247],
            [4, 1, 'M', 'HC', '_H3_', 0, 2.5, 0.03, 0.06, 1.425, 1.25,
             0.00859824, 0.268726247],
            [5, 1, 'M', 'HC', '_H4_', 0, 2.5, 0.03, 0.06, 1.425, 1.25,
             0.00859824, 0.268726247]]

        expected_bonds = [
            [1, 2, 340.0, 1.09],
            [1, 3, 340.0, 1.09],
            [1, 4, 340.0, 1.09],
            [1, 5, 340.0, 1.09]]

        expected_angles = [
            [2, 1, 3, 33.0, 107.8],
            [2, 1, 4, 33.0, 107.8],
            [2, 1, 5, 33.0, 107.8],
            [3, 1, 4, 33.0, 107.8],
            [3, 1, 5, 33.0, 107.8],
            [4, 1, 5, 33.0, 107.8]]

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_nonbonding,
                         expected_bonds=expected_bonds,
                         expected_angles=expected_angles)

        # Load molecule
        opls2005 = OPLS2005ForceField()
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')

        # Workaround to parameterize with OPLS2005 without the Schrodinger
        # dependency
        parameters = parameterize_opls2005(opls2005, molecule, ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_propers = [
            [3, 1, 2, 5, 7.0, -1, 2, 0.0],
            [3, 1, 2, 6, 7.0, -1, 2, 0.0],
            [4, 1, 2, 5, 7.0, -1, 2, 0.0],
            [4, 1, 2, 6, 7.0, -1, 2, 0.0]]

        expected_impropers = [
            [3, 4, 1, 2, 15.0, -1, 2],
            [5, 6, 2, 1, 15.0, -1, 2]]

        # Check it up
        check_parameters(topology,
                         expected_propers=expected_propers,
                         expected_impropers=expected_impropers)

    def test_openffopls2005_parameterizer(self):
        """It checks the parameterized method."""

        from peleffy.topology import Molecule
        from peleffy.topology import Topology
        from peleffy.forcefield import OpenFFOPLS2005ForceField
        from peleffy.utils import get_data_file_path
        from .utils import check_parameters, parameterize_openffopls2005

        FORCE_FIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

        # Load molecule
        hybridff = OpenFFOPLS2005ForceField(FORCE_FIELD_NAME)
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        ffld_file = get_data_file_path('tests/MET_ffld_output.txt')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_off_nonbonding = [
            [1, 0, 'M', 'OFFT', '_C1_', 0, 3.3996695084235347, 0.1094,
             -0.1088, 0, 1.6998347542117673, 0, 0],
            [2, 1, 'M', 'OFFT', '_H1_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [3, 1, 'M', 'OFFT', '_H2_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [4, 1, 'M', 'OFFT', '_H3_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [5, 1, 'M', 'OFFT', '_H4_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0]]

        expected_off_bonds = [
            [1, 2, 376.8940758588, 1.094223427522],
            [1, 3, 376.8940758588, 1.094223427522],
            [1, 4, 376.8940758588, 1.094223427522],
            [1, 5, 376.8940758588, 1.094223427522]]

        expected_off_angles = [
            [2, 1, 3, 33.78875634641, 110.2468561538],
            [2, 1, 4, 33.78875634641, 110.2468561538],
            [2, 1, 5, 33.78875634641, 110.2468561538],
            [3, 1, 4, 33.78875634641, 110.2468561538],
            [3, 1, 5, 33.78875634641, 110.2468561538],
            [4, 1, 5, 33.78875634641, 110.2468561538]]

        expected_opls_nonbonding = [
            [1, 0, 'M', 'CT', '_C1_', 0, 3.5, 0.066, -0.1088, 1.975, 1.75,
             0.005, -0.74168571],
            [2, 1, 'M', 'HC', '_H1_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [3, 1, 'M', 'HC', '_H2_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [4, 1, 'M', 'HC', '_H3_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [5, 1, 'M', 'HC', '_H4_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247]]

        expected_opls_bonds = [
            [1, 2, 340.0, 1.09],
            [1, 3, 340.0, 1.09],
            [1, 4, 340.0, 1.09],
            [1, 5, 340.0, 1.09]]

        expected_opls_angles = [
            [2, 1, 3, 33.0, 107.8],
            [2, 1, 4, 33.0, 107.8],
            [2, 1, 5, 33.0, 107.8],
            [3, 1, 4, 33.0, 107.8],
            [3, 1, 5, 33.0, 107.8],
            [4, 1, 5, 33.0, 107.8]]

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_off_nonbonding,
                         expected_bonds=expected_off_bonds,
                         expected_angles=expected_off_angles)

        # Hybrid force field set up
        hybridff.set_nonbonding_parameters('OPLS2005')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_opls_nonbonding,
                         expected_bonds=expected_off_bonds,
                         expected_angles=expected_off_angles)

        # Hybrid force field set up
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_bond_parameters('OPLS2005')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_off_nonbonding,
                         expected_bonds=expected_opls_bonds,
                         expected_angles=expected_off_angles)

        # Hybrid force field set up
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_bond_parameters('OpenFF')
        hybridff.set_angle_parameters('OPLS2005')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_off_nonbonding,
                         expected_bonds=expected_off_bonds,
                         expected_angles=expected_opls_angles)

        # Hybrid force field set up
        hybridff.set_nonbonding_parameters('OPLS2005')
        hybridff.set_bond_parameters('OPLS2005')
        hybridff.set_angle_parameters('OPLS2005')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_opls_nonbonding,
                         expected_bonds=expected_opls_bonds,
                         expected_angles=expected_opls_angles)

        # Load molecule
        hybridff = OpenFFOPLS2005ForceField(FORCE_FIELD_NAME)
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Define expected parameters
        expected_off_nonbonding = [
            [1, 0, 'M', 'OFFT', '_C1_', 0, 3.3996695084235347, 0.086,
             -0.218, 0, 1.6998347542117673, 0, 0],
            [2, 1, 'M', 'OFFT', '_C2_', 0, 3.3996695084235347, 0.086,
             -0.218, 0, 1.6998347542117673, 0, 0],
            [3, 1, 'M', 'OFFT', '_H1_', 0, 2.59964245953351, 0.015,
             0.109, 0, 1.299821229766755, 0, 0],
            [4, 1, 'M', 'OFFT', '_H2_', 0, 2.59964245953351, 0.015,
             0.109, 0, 1.299821229766755, 0, 0],
            [5, 2, 'M', 'OFFT', '_H3_', 0, 2.59964245953351, 0.015,
             0.109, 0, 1.299821229766755, 0, 0],
            [6, 2, 'M', 'OFFT', '_H4_', 0, 2.59964245953351, 0.015,
             0.109, 0, 1.299821229766755, 0, 0]]

        expected_off_bonds = [[1, 2, 404.83941263565, 1.372230219368],
                              [1, 3, 404.20804685, 1.085503378387],
                              [1, 4, 404.20804685, 1.085503378387],
                              [2, 5, 404.20804685, 1.085503378387],
                              [2, 6, 404.20804685, 1.085503378387]]

        expected_off_angles = [[1, 2, 5, 34.202963712735, 133.1339832262],
                               [1, 2, 6, 34.202963712735, 133.1339832262],
                               [2, 1, 3, 34.202963712735, 133.1339832262],
                               [2, 1, 4, 34.202963712735, 133.1339832262],
                               [3, 1, 4, 27.86109944498, 134.0642036482],
                               [5, 2, 6, 27.86109944498, 134.0642036482]]

        expected_off_propers = [
            [3, 1, 2, 5, 5.376019778605, -1, 2, 0.0],
            [3, 1, 2, 6, 5.376019778605, -1, 2, 0.0],
            [4, 1, 2, 5, 5.376019778605, -1, 2, 0.0],
            [4, 1, 2, 6, 5.376019778605, -1, 2, 0.0]]

        expected_off_impropers = [
            [1, 2, 5, 6, 1.1, -1, 2],
            [2, 1, 3, 4, 1.1, -1, 2]]

        expected_opls_propers = [
            [3, 1, 2, 5, 7.0, -1, 2, 0.0],
            [3, 1, 2, 6, 7.0, -1, 2, 0.0],
            [4, 1, 2, 5, 7.0, -1, 2, 0.0],
            [4, 1, 2, 6, 7.0, -1, 2, 0.0]]

        expected_opls_impropers = [
            [3, 4, 1, 2, 15.0, -1, 2],
            [5, 6, 2, 1, 15.0, -1, 2]]

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_off_nonbonding,
                         expected_bonds=expected_off_bonds,
                         expected_angles=expected_off_angles,
                         expected_propers=expected_off_propers,
                         expected_impropers=expected_off_impropers)

        # Hybrid force field set up
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_dihedral_parameters('OPLS2005')

        # Workaround to skip Schrodinger dependency
        parameters = parameterize_openffopls2005(hybridff, molecule,
                                                 ffld_file)

        # Generate molecular topology
        topology = Topology(molecule, parameters)

        # Check it up
        check_parameters(topology,
                         expected_nonbonding=expected_off_nonbonding,
                         expected_bonds=expected_off_bonds,
                         expected_angles=expected_off_angles,
                         expected_propers=expected_opls_propers,
                         expected_impropers=expected_opls_impropers)
