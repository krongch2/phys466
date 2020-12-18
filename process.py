import os
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import seaborn as sns

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
colors = sns.color_palette()

thz_to_cm1 = 33.35641

def post():

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
    # phonon.plot_band_structure().savefig('band.png')

    def myband():
        d = phonon.get_band_structure_dict()
        q = d['qpoints']
        freq = d['frequencies']
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
        x = q[-1]
        y = freq[-1]
        for i in range(y.shape[1]):
            ax.plot(x, y[:, i]*thz_to_cm1, color=colors[0])
        xlim = [0, 1/3]
        ax.set(
            xticklabels=labels, xticks=xlim, xlim=xlim,
            ylabel='Phonon frequency~$\\mathrm{cm}^{-1}$'
        )
        fig.tight_layout()
        plt.savefig('myband.png', bbox_inches='tight')
        os.system('rsub myband.png')

    myband()
    freq, eigvecs = phonon.get_frequencies_with_eigenvectors(path[0])
    print(freq)
    print(eigvecs)

    print(phonon.get_mesh_dict())

    phonon.set_total_DOS()
    phonon.get_total_dos_dict
    phonon.plot_total_DOS().savefig('dos.png')

    phonon.set_thermal_properties()
    tp_dict = phonon.get_thermal_properties_dict()
    with open('tp.pkl', 'wb') as f:
        pickle.dump(tp_dict, f)
    phonon.plot_thermal_properties().savefig('therm.png')

post()
