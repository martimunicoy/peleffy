"""
This module contains tests that check that the BCEConformations module.
"""

import numpy as np
import pytest


def build_mock_BCEConformations(pdb_file, ffld_file):
    from peleffy.topology import Molecule, BCEConformations
    from peleffy.forcefield import ForceFieldSelector
    from peleffy.topology import Topology
    from .utils import parameterize_opls2005

    molecule = Molecule(pdb_file)
    ff_selector = ForceFieldSelector()
    forcefield = ff_selector.get_by_name("opls2005")
    parameters = parameterize_opls2005(forcefield, molecule, ffld_file, charge_method="opls2005")
    topology = Topology(molecule, parameters)
    return BCEConformations(topology, "")


class TestBCEConformations(object):
    """BCEConformations test."""

    from peleffy.topology import BCEConformations

    def test_calculate_cluster_offsets(self):
        from peleffy.utils import get_data_file_path

        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        ffld_path = get_data_file_path("tests/ETL_ffld_output.txt")
        golden_offsets = [["_C1_", 0.0, 0.0, 0.0],
                          ["_C2_", -1.305,  0.181, -0.014],
                          ["_H1_", 0.72, 0.822, 0.042],
                          ["_H2_", 0.448, -0.988, -0.029],
                          ["_H3_", -1.73, 1.171, 0.015],
                          ["_H4_", -1.936, -0.687, -0.056]]
        bce_obj = build_mock_BCEConformations(pdb_path, ffld_path)
        bce_obj.calculate_cluster_offsets(pdb_path)
        for dih1, dih2 in zip(bce_obj.conformations_library[pdb_path], golden_offsets):
            assert dih1[0] == dih2[0]
            np.testing.assert_almost_equal(dih1[1:], dih2[1:], decimal=3)

    def test_write_dihedral_library(self):
        import os
        import tempfile
        from peleffy.utils import get_data_file_path, temporary_cd

        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        ffld_path = get_data_file_path("tests/ETL_ffld_output.txt")
        bce_obj = build_mock_BCEConformations(pdb_path, ffld_path)
        bce_obj.calculate_cluster_offsets(pdb_path)
        calculated_lines = []
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                bce_obj.save(os.path.join(tmpdir, "ETH.conformation"))
                with open(os.path.join(tmpdir, "ETH.conformation")) as f:
                    calculated_lines = f.readlines()
        calculated_lines = calculated_lines[3:]
        golden_conformation_library_path = get_data_file_path('parameters/ETH.conformation')
        with open(golden_conformation_library_path) as f:
            golden_lines = f.readlines()[3:]

        assert golden_lines == calculated_lines
