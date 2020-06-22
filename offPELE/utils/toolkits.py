import importlib
from distutils.spawn import find_executable
import tempfile
import os
import subprocess

import numpy as np
from simtk import unit

from offPELE.utils import temporary_cd


# Exceptions
class ToolkitUnavailableException(Exception):
    """The requested toolkit is unavailable."""
    pass


class ChargeCalculationError(Exception):
    """An external error when calculating charges"""
    pass


class ToolkitWrapper(object):
    """
    Toolkit wrapper base class.
    """
    _is_available = None
    _toolkit_name = None

    @property
    def toolkit_name(self):
        return self._toolkit_name

    @staticmethod
    def is_available():
        """
        Check whether the corresponding toolkit can be imported
        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.
        """
        return NotImplementedError


class RDKitToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'RDKit Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

    @staticmethod
    def is_available():
        """
        Check whether the RDKit toolkit can be imported
        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.
        """
        try:
            importlib.import_module('rdkit', 'Chem')
            return True
        except ImportError:
            return False

    def to_sfd_file(self, molecule, path):
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule

        with open('molecule.sdf', 'w') as f:
            writer = Chem.SDWriter(f)
            writer.write(rdkit_molecule)
            writer.close()


class AmberToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'Amber Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

        self._rdkit_toolkit_wrapper = RDKitToolkitWrapper()

    @staticmethod
    def is_available():
        """
        Check whether the AmberTools toolkit is installed
        Returns
        -------
        is_installed : bool
            True if AmberTools is installed, False otherwise.
        """
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            return False
        if not(RDKitToolkitWrapper.is_available()):
            return False
        return True

    def compute_partial_charges_am1bcc(self, molecule):
        off_molecule = molecule.off_molecule

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = off_molecule.total_charge / \
                    unit.elementary_charge

                self._rdkit_toolkit_wrapper.to_sfd_file(molecule, tmpdir)

                subprocess.check_output([
                    "antechamber", "-i", "molecule.sdf", "-fi", "sdf",
                    "-o", "charged.ac", "-fo", "ac", "-pf", "yes", "-dr", "n",
                    "-c", "bcc", "-nc", str(net_charge)])
                # Write out just charges
                subprocess.check_output([
                    "antechamber", "-dr", "n", "-i", "charged.ac", "-fi", "ac",
                    "-o", "charged2.ac", "-fo", "ac", "-c", "wc",
                    "-cf", "charges.txt", "-pf", "yes"])

                if not os.path.exists('charges.txt'):
                    # TODO: copy files into local directory to aid debugging?
                    raise ChargeCalculationError(
                        "Antechamber/sqm partial charge calculation failed on "
                        "molecule {} (SMILES {})".format(
                            off_molecule.name, off_molecule.to_smiles()))

                # Read the charges
                with open('charges.txt', 'r') as infile:
                    contents = infile.read()

                text_charges = contents.split()
                charges = np.zeros([molecule.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)

        charges = unit.Quantity(charges, unit.elementary_charge)


class OpenForceFieldToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'OpenForceField Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

    @staticmethod
    def is_available():
        """
        Check whether the OpenForceField toolkit is installed
        Returns
        -------
        is_installed : bool
            True if OpenForceField is installed, False otherwise.
        """
        try:
            importlib.import_module('openforcefield')
            return True
        except ImportError:
            return False
