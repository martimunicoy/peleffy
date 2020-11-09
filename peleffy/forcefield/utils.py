from peleffy.utils import get_data_file_path
from peleffy.topology.molecule import Molecule
from peleffy.utils.toolkits import RDKitToolkitWrapper

def checkBonds(scale, atom_type, degree, parentH):
    """
    """
    if atom_type == 'H' and parentH == 'O': scale = '1.05'
    if atom_type == 'H' and parentH == 'N': scale = '1.15'
    if atom_type == 'C' and degree == 3: scale= '1.875'
    if atom_type == 'C' and degree == 2: scale = '1.825'
    if atom_type == 'N' and degree == 4: scale = '1.625'
    if atom_type == 'N' and degree == 1: scale = '1.60'
    if atom_type == 'O' and degree == 1: scale = '1.48'
    return scale

def find_GBSA_parameters_according_to(atom_name, atom_type, degree,
                                      parentH):
    # Parameter List extracted from PELE solvent templates for OBC
    paramtersLst = [['CW', '1.875', '0.72'],
                        ['NC', '1.7063', '0.79'],
                        ['CM', '1.875', '0.72'],
                        ['C*', '1.875', '0.72'],
                        ['H1', '1.25', '0.85'],
                        ['CT', '1.9', '0.72'],
                        ['N2', '1.7063', '0.79'],
                        ['N*', '1.7063', '0.79'],
                        ['CR', '1.875', '0.72'],
                        ['HO', '1.05', '0.85'],
                        ['NB', '1.7063', '0.79'],
                        ['H2', '1.25', '0.85'],
                        ['S', '1.775', '0.96'],
                        ['NA', '1.7063', '0.79'],
                        ['H4', '1.25', '0.85'],
                        ['HC', '1.25', '0.85'],
                        ['C', '1.875', '0.72'],
                        ['OH', '1.535', '0.85'],
                        ['CQ', '1.875', '0.72'],
                        ['CK', '1.875', '0.72'],
                        ['O2', '1.48', '0.85'],
                        ['OS', '1.535', '0.85'],
                        ['SH', '1.775', '0.96'],
                        ['HA', '1.25', '0.85'],
                        ['CB', '1.875', '0.72'],
                        ['H5', '1.25', '0.85'],
                        ['CN', '1.875', '0.72'],
                        ['P', '1.87', '0.86'],
                        ['N3', '1.625', '0.79'],
                        ['HP', '1.25', '0.85'],
                        ['N', '1.7063', '0.79'],
                        ['H', '1.15', '0.85'],
                        ['HS', '1.25', '0.85'],
                        ['CV', '1.875', '0.72'],
                        ['CA', '1.875', '0.72'],
                        ['O', '1.48', '0.85'],
                        ['CC', '1.875', '0.72'], ]
    #  ['HWS', '1.05', '0.85'] ,
    #  ['OWS', '1.535', '0.85'] ,]
    # Parameter List extracted from TINKER source code for OBC
    # (ksolv.f: +/-line 423)
    atomTypesOverlapFactors = [['H', '1.25'],
                               ['Li', '1.432'],
                               ['C', '1.90'],
                               ['N', '1.7063'],
                               ['O', '1.535'],
                               ['F', '1.47'],
                               ['FE', '2.00'],  # default parameters
                               ['Ne', '1.39'],
                               ['Na', '1.992'],
                               ['Mg', '1.70'],
                               ['Si', '1.80'],
                               ['P', '1.87'],
                               ['S', '1.775'],
                               ['Cl', '1.735'],
                               ['Ar', '1.70'],
                               ['K', '2.123'],
                               ['Ca', '1.817'],
                               ['Br', '1.90'],
                               ['Kr', '1.812'],
                               ['Rb', '2.26'],
                               ['I', '2.10'],
                               ['Xe', '1.967'],
                               ['Cs', '2.507'],
                               ['Ba', '2.188'],
                               ['Pt', '2.0'],
                               ]
    # Parameter List extracted from TINKER source code for OBC
    # (ksolv.f: +/-line 423)
    atomTypesHCTradii = [['H', '0.85'],
                         ['C', '0.72'],
                         ['N', '0.79'],
                         ['O', '0.85'],
                         ['F', '0.88'],
                         ['P', '0.86'],
                         ['S', '0.96'],
                         ['Pt', '0.80'],
                         ['FE', '0.88'],
                         ]

    # Assign overlapFactors and HCT radii using the atom name
    for param in paramtersLst:
        if atom_name.strip() == param[0]:
            scale, radius= param[1], param[2]
            break

    # Assign overlapFactors and HCT radii using the atom type
    else:
        found = False
        for OverlapFactor in atomTypesOverlapFactors:
            if atom_type.strip() == OverlapFactor[0]:
                scale, found = OverlapFactor[1], True
                scale = checkBonds(scale, atom_type, degree, parentH)
                break
        for HCTradii in atomTypesHCTradii:
            if atom_type.strip() == HCTradii[0]:
                radius, found = HCTradii[1], (found and True)
                break
        else: found = False

    # Retuns overlapFractor and HCT radii if found otherwise it raises a
    # warining and returns the default parameters
        if not found:
            from peleffy.utils import Logger
            log = Logger()
            log.warning('Parameter for {} {} NOT found in the template '.format(
                  atom_name,atom_type) + 'database... using default parameters')
            radius, scale = '0.80','2.0'
    return radius, scale




