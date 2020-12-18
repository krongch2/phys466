import os
import pickle

import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

def run_phonon():

    phonon = phonopy.load(
        supercell_matrix=[2, 2, 2],
        primitive_matrix='auto',
        unitcell_filename="POSCAR_unitcell",
        force_constants_filename="FORCE_CONSTANTS"
    )

    phonon.set_mesh([20, 20, 20])
    path = [[0, 0, 0], [2/3, 1/3, 0]]
    labels = ["$\\Gamma$", "K"]
    qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
    phonon.run_band_structure(qpoints, path_connections=connections, labels=labels, with_eigenvectors=True)
    phonon.set_total_DOS()
    phonon.set_thermal_properties()
    pp = {
        'band': phonon.get_band_structure_dict(),
        'freq_eigvecs': phonon.get_frequencies_with_eigenvectors(path[0]),
        'dos': phonon.get_total_dos_dict(),
        'therm': phonon.get_thermal_properties_dict()
    }
    with open('pp.pkl', 'wb') as f:
        pickle.dump(pp, f)

run_phonon()
