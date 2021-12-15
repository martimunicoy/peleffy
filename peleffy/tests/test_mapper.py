"""
This module contains the tests to check peleffy's molecular mapper.
"""

import pytest


class TestMapper(object):
    """
    It wraps all tests that involve Mapper class.
    """

    def test_mapper_initializer(self):
        """
        It checks the initialization of the Mapper class.
        """
        from peleffy.topology import Molecule
        from peleffy.topology import Mapper

        mol1 = Molecule(smiles='c1ccccc1', hydrogens_are_explicit=False)
        mol2 = Molecule(smiles='c1ccccc1C', hydrogens_are_explicit=False)

        # Check initializer with only the two molecules
        mapper = Mapper(mol1, mol2)

        # Check initializer with only include_hydrogens parameter
        mapper = Mapper(mol1, mol2, include_hydrogens=False)

        # Check initializer with bad types
        with pytest.raises(TypeError):
            mapper = Mapper(mol1.rdkit_molecule, mol2)

        with pytest.raises(TypeError):
            mapper = Mapper(mol1, "mol2")

    def test_mapper_mapping(self):
        """
        It validates the mapping.
        """
        from peleffy.topology import Molecule
        from peleffy.topology import Mapper

        # First mapping checker
        mol1 = Molecule(smiles='c1ccccc1', hydrogens_are_explicit=False)
        mol2 = Molecule(smiles='c1ccccc1C', hydrogens_are_explicit=False)

        mapper = Mapper(mol1, mol2, include_hydrogens=False)
        mapping = mapper.get_mapping()

        assert mapping == [(0, 0), (1, 1), (2, 2), (3, 3),
                           (4, 4), (5, 5)], 'Unexpected mapping'

        # Second mapping checker
        mol1 = Molecule(smiles='c1(C)ccccc1C', hydrogens_are_explicit=False)
        mol2 = Molecule(smiles='c1c(C)cccc1C', hydrogens_are_explicit=False)

        mapper = Mapper(mol1, mol2, include_hydrogens=False)
        mapping = mapper.get_mapping()

        assert mapping == [(0, 1), (1, 2), (2, 0), (3, 6),
                           (4, 5), (5, 4), (6, 3)], 'Unexpected mapping'

        # Third mapping checker with hydrogens
        mol1 = Molecule(smiles='c1ccccc1', hydrogens_are_explicit=False)
        mol2 = Molecule(smiles='c1ccccc1C', hydrogens_are_explicit=False)

        mapper = Mapper(mol1, mol2, include_hydrogens=True)
        mapping = mapper.get_mapping()

        assert mapping == [(0, 0), (1, 1), (2, 2), (3, 3),
                           (4, 4), (5, 5), (11, 6), (10, 11),
                           (9, 10), (8, 9), (7, 8), (6, 7)], \
            'Unexpected mapping'

        # Fourth mapping checker with hydrogens
        mol1 = Molecule(smiles='c1(C)ccccc1C', hydrogens_are_explicit=False)
        mol2 = Molecule(smiles='c1c(C)cccc1C', hydrogens_are_explicit=False)

        mapper = Mapper(mol1, mol2, include_hydrogens=True)
        mapping = mapper.get_mapping()

        assert mapping == [(0, 1), (1, 2), (8, 9), (9, 10),
                           (10, 11), (2, 0), (3, 6), (4, 5),
                           (5, 4), (6, 3), (7, 12), (14, 13),
                           (13, 14), (12, 7), (11, 8)], 'Unexpected mapping'
