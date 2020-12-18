import os
import pickle
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
colors = sns.color_palette()

thz_to_cm1 = 33.35641

class tblg_lammps_relaxed:
    """attempting to make data extraction cleaner"""
    def __init__(self, filename):
        filename_split = filename.split('/')
        # angle = filename[-2]
        # angle = float('0.'+angle[5:])

        data_raw = open(filename, 'r')
        data_lines = data_raw.readlines()
        N = int(data_lines[2].split(' ')[0])
        self.N = N
        #now we're just extracting the data....
        xlo,xhi = data_lines[5].split(' ')[0], data_lines[5].split(' ')[1]
        self.x_axis = float(xhi) - float(xlo)
        ylo,yhi = data_lines[6].split(' ')[0], data_lines[6].split(' ')[1]
        self.y_axis = float(yhi) - float(ylo)
        zlo,zhi = data_lines[7].split(' ')[0], data_lines[7].split(' ')[1]
        self.z_axis = float(zhi) - float(zlo)
        self.xy_tilt = float(data_lines[8].split(' ')[0])

        #atomic coordinates...
        self.coords = np.zeros([N,3])
        for i in range(N):
            xi,yi,zi = data_lines[16+i].split(' ')[2:5]
            self.coords[i,0] = float(xi)
            self.coords[i,1] = float(yi)
            self.coords[i,2] = float(zi)
        halfN = int(N/2)
        self.base_height = np.mean(self.coords[0:halfN,2])
        self.base_std = np.std(self.coords[0:halfN,2])
        self.twist_height = np.mean(self.coords[halfN:,2])
        self.twist_std = np.std(self.coords[halfN:,2])
        self.dz = self.twist_height -self.base_height
        self.base_buckled = 0
        self.twist_buckled = 0
        for i in range(halfN):
            # if self.coords[i,2] != self.coords[0,2]: self.base_buckled +=1
            if abs(self.coords[i,2] - self.base_height)> 1e-6: self.base_buckled +=1
            # if self.coords[i+halfN,2] != self.coords[0,2]: self.twist_buckled +=1
            if abs(self.coords[i+halfN,2] - self.twist_height)> 1e-6: self.twist_buckled +=1
        return

def get_natoms():
    with open('POSCAR_unitcell') as f:
        text = f.read()
        for line in text.split('\n'):
            m =  re.search('^\d+$', line)
            if m:
                natoms = int(m.group(0))
    return natoms

def plot_band(band, pot, angle):
    q = band['qpoints']
    freq = band['frequencies']
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    x = q[-1]
    y = freq[-1]
    for i in range(y.shape[1]):
        ax.plot(x, y[:, i]*thz_to_cm1, color=colors[0])
    xlim = [0, 1/3]
    ax.set(
        xticklabels=['$\\Gamma$', 'K'], xticks=xlim, xlim=xlim,
        ylabel='Phonon frequency~$\\mathrm{cm}^{-1}$',
        title=f'{pot} {angle} deg'
    )
    fig.tight_layout()
    plt.savefig('band.png', bbox_inches='tight')
    os.system('rsub band.png')

# def plot_dos(dos):
#     fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
#     ax.plot(dos['frequency_points'], dos['total_dos']/get_natoms(), '-', color=colors[0])
#     fig.tight_layout()
#     plt.savefig('dos.png', bbox_inches='tight')
#     os.system('rsub dos.png')

def plot_therm(therm):
    # entropy, free_energy
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    ax.plot(therm['temperatures'], therm['heat_capacity']/get_natoms(), '-', color=colors[0])
    ax.set(
        xlabel='$T$',
        ylabel='$C$')
    fig.tight_layout()
    plt.savefig('heat.png', bbox_inches='tight')
    os.system('rsub heat.png')

def plot_eigvecs(freq_eigvecs):
    freq, eigvecs = freq_eigvecs
    print(freq)
    print(eigvecs)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    # pos = tblg_lammps_relaxed()
    # ax.plot()

def read_pp():

    heat_l = []
    dos_l = []
    workdir = os.getcwd()
    os.chdir('data')
    for dirname in sorted(os.listdir('.')):
        os.chdir(dirname)
        # and (dirname == 'airebo_13.2' or dirname == 'airebo_21.7')
        if os.path.isfile('pp.pkl'):
            print(dirname)
            pot, angle = dirname.split('_')
            with open('pp.pkl', 'rb') as f:
                pp = pickle.load(f)
            plot_band(pp['band'], pot, angle)

            # plot_eigvecs(pp['freq_eigvecs'])
            dos = pd.DataFrame(pp['dos'])
            dos['pot'] = pot
            dos['angle'] = angle
            dos['natoms'] = get_natoms()
            dos_l.append(dos)

            heat = pd.DataFrame(pp['therm'])
            heat['pot'] = pot
            heat['angle'] = angle
            heat['natoms'] = get_natoms()
            heat_l.append(heat)

        os.chdir('..')
    os.chdir('..')
    pd.concat(heat_l, ignore_index=True).to_csv('heat.csv', index=False)
    pd.concat(dos_l, ignore_index=True).to_csv('dos.csv', index=False)

def plot_heat():
    d = pd.read_csv('heat.csv')
    d['C_per_atom'] = d['heat_capacity'] / d['natoms']
    d['S_per_atom'] = d['entropy'] / d['natoms']
    d['F_per_atom'] = d['free_energy'] / d['natoms']

    def subtract(g):
        mask = g['angle'] == 0.0
        g['dC'] = g['C_per_atom'] - float(g.loc[mask, 'C_per_atom'])
        g['dS'] = g['S_per_atom'] - float(g.loc[mask, 'S_per_atom'])
        g['dF'] = g['F_per_atom'] - float(g.loc[mask, 'F_per_atom'])
        return g

    p = d.groupby(['temperatures', 'pot']).apply(subtract)
    print(p)

    g = sns.FacetGrid(p, col='pot', hue='angle', height=3, sharey=False)
    g.map(sns.scatterplot, 'temperatures', 'dC')
    for ax in g.axes.flatten():
        ax.set(
            xlabel='$T~(\\mathrm{K})$',
            ylabel='$\\Delta C~(\\mathrm{J/K})$'
        )

    g.add_legend()
    output = 'heat.png'
    plt.savefig(f'{output}', bbox_inches='tight')
    os.system(f'rsub {output}')

def plot_dos():
    d = pd.read_csv('dos.csv')
    d['dos'] = d['total_dos'] / d['natoms']
    d['omega'] = d['frequency_points']*thz_to_cm1

    g = sns.FacetGrid(d, col='pot', hue='angle', height=3, sharex=False, sharey=False)
    g.map(sns.lineplot, 'omega', 'dos')
    for ax in g.axes.flatten():
        ax.set(
            xlabel='$\\omega~(\\mathrm{cm}^{-1})$',
            ylabel='$D(\\omega)$'
        )
    g.add_legend()
    output = 'dos.png'
    plt.savefig(f'{output}', bbox_inches='tight')
    os.system(f'rsub {output}')

read_pp()
# plot_heat()
# plot_dos()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All you need to do for this is call:

variable_for_relaxed_angle_A("path-to-relax.pos")

this only reads the output because the line spacing on the input is slightly different
@author: gnolan2
"""
# pos = tblg_lammps_relaxed('data/airebo_21.7/relax.pos')

# print(pos.coords)
# a061_relax = tblg_lammps_relaxed('angle061/data.relaxed_CG')
