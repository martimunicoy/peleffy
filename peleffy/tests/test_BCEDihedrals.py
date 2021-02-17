"""
This module contains tests that check that the BCEDihedrals module.
"""

import numpy as np
import pytest

def build_mock_BCEDihedrals(pdb_file):
    from peleffy.topology import Molecule
    from peleffy.BCEDihedrals import BCEDihedrals
    from peleffy.forcefield import ForceFieldSelector
    from peleffy.topology import Topology

    molecule = Molecule(pdb_file)
    ff_selector = ForceFieldSelector()
    forcefield = ff_selector.get_by_name("opls2005")
    parameters = forcefield.parameterize(molecule, charge_method="opls2005")
    topology = Topology(molecule, parameters)
    return BCEDihedrals(topology, "", mode="")

class TestBCEDihderals(object):
    """BCEDihedrals test."""

    from peleffy.BCEDihedrals import BCEDihedrals

    def test_list_dihedrals(self):
        from peleffy.utils import get_data_file_path

        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        golden_dihedrals = [(2, 0, 1, 4), (2, 0, 1, 5), (3, 0, 1, 4), (3, 0, 1, 5)]
        bce_obj = build_mock_BCEDihedrals(pdb_path)
        dihedrals_indices = bce_obj.list_all_dihedrals()
        dihedrals_indices.sort()
        np.testing.assert_array_equal(dihedrals_indices, golden_dihedrals)

    def test_calculate_cluster_angles(self):
        from peleffy.utils import get_data_file_path

        golden_dihedrals = [(2, 0, 1, 4), (2, 0, 1, 5), (3, 0, 1, 4), (3, 0, 1, 5)]
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        golden_angles = [ ["_H1_", "_C1_", "_C2_", "_H3_", 0.001151],
                          ["_H1_", "_C1_", "_C2_", "_H4_", -3.141278],
                          ["_H2_", "_C1_", "_C2_", "_H3_", -3.141278],
                          ["_H2_", "_C1_", "_C2_", "_H4_", -0.000366]]
        bce_obj = build_mock_BCEDihedrals(pdb_path)
        bce_obj.calculate_cluster_angles(pdb_path, golden_dihedrals, match_indexes=False)
        for dih1, dih2 in zip(bce_obj.dihedral_library[pdb_path], golden_angles):
            np.testing.assert_array_equal(dih1[:4], dih2[:4])
            np.testing.assert_almost_equal(dih1[-1], dih2[-1], decimal=3)

    def test_write_dihedral_library(self):
        import os
        import tempfile
        from peleffy.utils import get_data_file_path, temporary_cd

        golden_dihedrals = [(2, 0, 1, 4), (2, 0, 1, 5), (3, 0, 1, 4), (3, 0, 1, 5)]
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        bce_obj = build_mock_BCEDihedrals(pdb_path)
        bce_obj.calculate_cluster_angles(pdb_path, golden_dihedrals, match_indexes=False)
        calculated_lines = []
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                bce_obj.save(os.path.join(tmpdir, "ETH.dihedral"))
                with open(os.path.join(tmpdir, "ETH.dihedral")) as f:
                    calculated_lines = f.readlines()
        calculated_lines = calculated_lines[3:]
        golden_dihedral_library_path = get_data_file_path('parameters/ETH.dihedral')
        with open(golden_dihedral_library_path) as f:
            golden_lines = f.readlines()[3:]

        assert golden_lines == calculated_lines
