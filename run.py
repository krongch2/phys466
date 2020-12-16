import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from phonolammps import Phonolammps
import phonopy

def test():

    os.chdir('data')

    with open('lammps2.lsf', 'w') as f:
        f.write(
'''#!/bin/bash
#BSUB -P mat221
#BSUB -W 15
#BSUB -nnodes 1
#BSUB -J lmp
#BSUB -o lsf.%J
#BSUB -e lsf.%J

module load lammps/master

# jsrun --smpiargs="-gpu" -n1 -g6 -c42 -a6 --bind=proportional-packed:7 lmp -in data/dipole.in > data/dipole.out
jsrun --smpiargs="-gpu" -n1 -g6 -c42 -a6 --bind=proportional-packed:7 lmp -in dipole.in > dipole.out
# jsrun -n 1 -a 1 -g 1 -c 1 lmp -in dipole.in > dipole.out
''')

    os.system('bsub lammps2.lsf')
    os.chdir('..')


def writeLammpsData(f_coor,coor,a1,a2):
    f = open(f_coor,'w+')

    natoms = len(coor)

    xlim = [0,a1[0]]
    ylim = [0,a2[1]]
    zlim = [-20,20]
    xy = a2[0]
    xz = 0
    yz = 0

    f.write('%s (written by ASE)\n\n'%(f_coor))
    f.write('%d\tatoms\n'%(natoms))
    f.write('%d\tatom types\n'%(1))
    f.write('%.9g\t%.9g\txlo xhi\n'%(xlim[0],xlim[1]))
    f.write('%.9g\t%.9g\tylo yhi\n'%(ylim[0],ylim[1]))
    f.write('%.9g\t%.9g\tzlo zhi\n'%(zlim[0],zlim[1]))
    f.write('%.9g\t%.9g\t%.9g\txy xz yz\n'%(xy, xz, yz))

    f.write('\n\n')
    f.write('Masses\n\n')

    f.write('%d\t12\n'%(1))
    f.write('\n')

    f.write('Atoms\n\n')
    iatom = 0
    for j in range(len(coor)):
        iatom+=1
        f.write('%d %d %d %d %1.3f %1.3f %1.3f %d %d %d\n'%(iatom,1,1,0,coor[j,0],coor[j,1],coor[j,2],0,0,0))

    f.close()


def run_ml():

    # command_line = 'lmp_emil -var f_coor ' + flcoor + ' -var flout_all ' + flout + ' -var fl_energy ' + flenergy + ' < MLfindenergy.in'

    # print(command_line)
    # os.system(command_line)
    # coor = np.genfromtxt(flout, skip_header=2, usecols=[1, 2, 3])
    # teng, ecoh = np.genfromtxt(flenergy, delimiter=', ')

    # coor, teng = monolayerenergy(coords, a, b)
    os.system('cp ml.in data')
    os.system('cp potentials/CH.rebo data')
    os.system('cp potentials/CH_taper.KC data')

    os.chdir('data')

    a1 = 2.45275
    coords = a1*np.array([[0, 0, 0], [0.5, 3**0.5/6, 0]])


    a = np.array([a1, 0, 0])
    b = np.array([a1/2, 3**0.5*a1/2, 0])

    flcoor = 'coor_mlenergy.input'
    writeLammpsData(flcoor, coords, a, b)
    flout = 'coor_mlenergy.out'
    flenergy = 'energy_mlenergy.out'


#     with open('lammps.lsf', 'w') as f:
#         f.write(
# f'''#!/bin/bash
# #BSUB -P mat221
# #BSUB -W 15
# #BSUB -nnodes 1
# #BSUB -J lmp
# #BSUB -o lsf.%J
# #BSUB -e lsf.%J

# module load lammps/master

# # jsrun --smpiargs="-gpu" -n1 -g6 -c42 -a6 --bind=proportional-packed:7 lmp -in data/dipole.in > data/dipole.out
# # jsrun --smpiargs="-gpu" -n1 -g6 -c42 -a6 --bind=proportional-packed:7 lmp -in dipole.in > dipole.out
# jsrun -n 1 -a 1 -g 1 -c 1 lmp -var f_coor {flcoor} -var flout_all {flout} -var fl_energy {flenergy} -in ml.in > ml.out
# ''')
#     os.system('bsub lammps.lsf')
    os.system(f'mpirun -np 4 lmp_mpi -var f_coor {flcoor} -var flout_all {flout} -var fl_energy {flenergy} -in ml.in > ml.out')
    os.chdir('..')

# run_ml()

def run_phonolammps():
    n = 2
    phlammps = Phonolammps('in.graphene', supercell_matrix=[[n, 0, 0], [0, n, 0], [0, 0, n]])
    unitcell = phlammps.get_unitcell()
    force_constants = phlammps.get_force_constants()
    supercell_matrix = phlammps.get_supercell_matrix()
    print('unitcell')
    print('*'*50)
    print(unitcell)
    print('force const')
    print('*'*50)
    print(force_constants)
    print('supercell')
    print('*'*50)
    print(supercell_matrix)

    from phonopy import Phonopy
    phonon = Phonopy(unitcell, supercell_matrix)
    print(phonon)
    print(dir(phonon))



    phonon.set_force_constants(force_constants)
    phonon.set_mesh([20, 20, 20])

    # phonon.write_yaml_band_structure('band.conf')
    # phonon.set_band_structure('band.conf')
    # phonon.plot_band_structure().savefig('band.png')

    # phonon.set_total_DOS()
    # phonon.plot_total_DOS().savefig('dos.png')

    phonon.set_thermal_properties()
    # print(phonon.get_thermal_properties_dict())
    phonon.plot_thermal_properties().savefig('therm.png')

# run_phonolammps()

def postprocess():
    # os.system('phonolammps in.graphene --dim 2 2 2 -c POSCAR_unitcell')

    phonon = phonopy.load(
        supercell_matrix=[2, 2, 2],
        primitive_matrix='auto',
        unitcell_filename="POSCAR_unitcell",
        force_constants_filename="FORCE_CONSTANTS"
    )

    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections


    phonon.set_mesh([20, 20, 20])

    path = [[0, 0, 0], [0.5, 0.5, 0]]
    labels = ["$\\Gamma$", "K"]
    qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
    phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
    phonon.plot_band_structure()# .savefig('band.png')
    print(dir(phonon))

    def myband():
        d = phonon.get_band_structure_dict()
        q = d['qpoints']
        freq = d['frequencies']
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 3))
        print(len(q))
        x = q[-1]
        y = freq[-1]
        print(x)
        print(y)
        print('*'*50)
        for i in range(y.shape[1]):
            ax.plot(x, y[:, i], color='blue')


        plt.savefig('myband.png')
        os.system('rsub myband.png')

    myband()

    # phonon.auto_band_structure(plot=True).savefig('band.png')
    # print(phonon.get_mesh_dict())

    # phonon.set_total_DOS()
    # phonon.plot_total_DOS().show()

    phonon.set_thermal_properties()
    phonon.plot_thermal_properties().savefig('therm_cmd.png')


def preprocess():
    with open('tblg.slurm', 'w') as f:
        f.write(
'''#! /bin/bash
#SBATCH --job-name="tblg"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20

cd $SLURM_SUBMIT_DIR

phonolammps in.graphene --dim 2 2 2 -c POSCAR_unitcell
''')

    os.system('sbatch tblg.slurm')

# preprocess()
postprocess()
