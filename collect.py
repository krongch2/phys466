import os
import pickle
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import phonopy
import seaborn as sns

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

def read_fc():

    l = []
    os.chdir('data')
    for dirname in sorted(os.listdir('.')):

        os.chdir(dirname)
        os.system('rm tp.pkl')
        if os.path.isfile('FORCE_CONSTANTS') and not os.path.isfile('tp.pkl') and (dirname == 'airebo_00.0' or dirname == 'airebo_21.7'):
            print(dirname)
            with open('pp.slurm', 'w') as f:
                f.write(
f'''#! /bin/bash
#SBATCH --job-name="{dirname}"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20

cd $SLURM_SUBMIT_DIR

python ../../postprocess.py
''')
            # os.system('sbatch pp.slurm')
            os.system('python ../../postprocess.py')

        if os.path.isfile('tp.pkl'):
            with open('tp.pkl', 'rb') as f:
                tp_dict = pickle.load(f)

            with open('POSCAR_unitcell') as f:
                text = f.read()
                for line in text.split('\n'):
                    m =  re.search('^\d+$', line)
                    if m:
                        natoms = m.group(0)
            pot, angle = dirname.split('_')
            d = pd.DataFrame(tp_dict)
            d['pot'] = pot
            d['angle'] = float(angle)
            d['natoms'] = int(natoms)
            l.append(d)

        os.chdir('..')
    os.chdir('..')

    d = pd.concat(l, ignore_index=True)
    d.to_csv('tp.csv', index=False)

read_fc()

def plot(csv, output):
    d = pd.read_csv(csv)
    d['C_per_atom'] = d['heat_capacity'] / d['natoms']
    d['S_per_atom'] = d['entropy'] / d['natoms']
    d['F_per_atom'] = d['free_energy'] / d['natoms']

    def subtract(g):
        mask = g['angle'] == 4.4
        g['dC'] = g['C_per_atom'] - float(g.loc[mask, 'C_per_atom'])
        g['dS'] = g['S_per_atom'] - float(g.loc[mask, 'S_per_atom'])
        g['dF'] = g['F_per_atom'] - float(g.loc[mask, 'F_per_atom'])
        return g

    p = d.groupby(['temperatures', 'pot']).apply(subtract)
    print(p)

    g = sns.FacetGrid(p, row='pot', hue='angle', height=3)
    g.map(sns.scatterplot, 'temperatures', 'dF')
    for ax in g.axes.flatten():
        ax.set_xlim(0,200)
    g.add_legend()
    plt.savefig(f'{output}', bbox_inches='tight')
    os.system(f'rsub {output}')

# plot('tp.csv', 'tp.png')
