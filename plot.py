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

def mask(d, options):
    '''
    Selects data from the given options
    '''
    match = True
    for key in options.keys():
        if d[key].dtype == 'float64':
            d[key] = np.round(d[key], decimals=4)
            options[key] = np.round(options[key], decimals=4)
        match = match & (d[key] == options[key])
    return match

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

# def plot_band(band, pot, angle):
#     q = band['qpoints']
#     freq = band['frequencies']
#     fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
#     x = q[-1]
#     y = freq[-1]

#     for i in range(y.shape[1]):
#         ax.plot(x, y[:, i]*thz_to_cm1, color=colors[0])
#     xlim = [0, 1/3]
#     ax.set(
#         xticklabels=['$\\Gamma$', 'K'], xticks=xlim, xlim=xlim,
#         ylabel='Phonon frequency~$\\mathrm{cm}^{-1}$',
#         title=f'{pot} {angle} deg'
#     )
#     fig.tight_layout()
#     plt.savefig('band.png', bbox_inches='tight')
#     os.system('rsub band.png')

def collect_band(band):
    q = band['qpoints']
    freq = band['frequencies']
    x = q[-1]
    y = freq[-1]
    l = []
    for i in range(y.shape[1]):
        # print(len(x))
        # print(len(y[:, i]*thz_to_cm1))
        l.append(pd.DataFrame({
            'q': x,
            'freq': y[:, i]*thz_to_cm1,
            'band': i
        })),
    d = pd.concat(l, ignore_index=True)
    return d

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
    band_l = []
    workdir = os.getcwd()
    os.chdir('data')
    for dirname in sorted(os.listdir('.')):
        os.chdir(dirname)
        if os.path.isfile('pp.pkl'):
            print(dirname)
            pot, angle = dirname.split('_')
            with open('pp.pkl', 'rb') as f:
                pp = pickle.load(f)
            band = collect_band(pp['band'])
            band['pot'] = pot
            band['angle'] = angle
            print(band['angle'])
            band['natoms'] = get_natoms()
            band_l.append(band)

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
    pd.concat(band_l, ignore_index=True).to_csv('band.csv', index=False)

def plot_heat():
    d = pd.read_csv('heat.csv')
    d['C'] = d['heat_capacity'] / d['natoms']
    d['S'] = d['entropy'] / d['natoms']
    d['F'] = d['free_energy'] / d['natoms']

    def subtract(g):
        mask = g['angle'] == '00.0'
        g['dC'] = g['C'] - float(g.loc[mask, 'C'])
        g['dS'] = g['S'] - float(g.loc[mask, 'S'])
        g['dF'] = g['F'] - float(g.loc[mask, 'F'])
        return g

    p = d.groupby(['temperatures', 'pot']).apply(subtract)
    print(p)

    # g = sns.FacetGrid(p, col='pot', hue='angle', height=3, sharey=False)
    # g.map(sns.scatterplot, 'temperatures', 'dC')
    g = sns.FacetGrid(d, col='pot', hue='angle', height=3, sharey=False)
    g.map(sns.scatterplot, 'temperatures', 'C')
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
    d = d.loc[~d['angle'].isin(['aa', 'mo']), :]

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

def plot_band_mo():
    d = pd.read_csv('band.csv')
    options = {
        'angle': 'mo'
    }
    d['pot'] = d['pot'].map({'airebo': 'AIREBO', 'tersoff': 'Tersoff', 'lj': 'LJ'})
    # d = d.sort_values(by=['pot_sort'])


    mo = d.loc[mask(d, options), :]
    pots = ['AIREBO', 'Tersoff', 'LJ']
    fig, axes = plt.subplots(nrows=1, ncols=len(pots), figsize=(6, 3), sharey=True)
    for i, pot in enumerate(pots):
        ax = axes[i]
        dd = mo.loc[mo['pot'] == pot, :]
        for band in dd['band'].unique():
            ddd = dd.loc[dd['band'] == band, :]
            ax.plot(ddd['q'], ddd['freq'], '-', color=colors[0])
            xlim = [0, 1/3]
            ax.set(
                xticklabels=['$\\Gamma$', 'K'], xticks=xlim, xlim=xlim,
                title=f'{pot} Monolayer'
            )

    axes[0].set(ylabel='$\\omega~(\\mathrm{cm}^{-1})$')
    fig.tight_layout()
    output = 'band_mo.pdf'
    plt.savefig(f'{output}', bbox_inches='tight')
    if 'png' in output:
        os.system(f'rsub {output}')

def plot_band_ab_notwist():
    d = pd.read_csv('band.csv')
    options = {
        'angle': 0,
        'pot': 'airebo'
    }
    r = pd.read_csv('exp_ab.csv', header=None)
    r.columns = ['q', 'freq']
    print(r)
    dd = d.loc[mask(d, options), :]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    for band in dd['band'].unique():
        ddd = dd.loc[dd['band'] == band, :]
        ax.plot(ddd['q'], ddd['freq'], '-', color=colors[0])
        xlim = [0, 1/3]
        ax.set(
            xticklabels=['$\\Gamma$', 'K'], xticks=xlim, xlim=xlim,
            title=f'AIREBO, AB Stacking'
        )
    ax.plot(r['q'], r['freq'], 'x', color=colors[1])
    ax.set(ylabel='$\\omega~(\\mathrm{cm}^{-1})$')
    fig.tight_layout()
    output = 'band_ab.pdf'
    plt.savefig(f'{output}', bbox_inches='tight')
    if 'png' in output:
        os.system(f'rsub {output}')

def plot_band_ab():
    d = pd.read_csv('band.csv', dtype=str)
    options = {
        'pot': 'airebo'
    }
    r = pd.read_csv('exp_ab.csv', header=None)
    r.columns = ['q', 'freq']
    print(r)
    d = d.loc[mask(d, options), :]
    angles = d['angle'].unique()
    print(angles)
    # angles = ['07.3', '13.2', '21.7', 'aa']
    # print(angles)
    return
    fig, axes = plt.subplots(nrows=1, ncols=len(angles), figsize=(3, 3))
    for i, angle in enumerate(angles):
        ax = axes[i]
        dd = d.loc[d['angle'] == angle, :]
        for band in dd['band'].unique():
            ddd = dd.loc[dd['band'] == band, :]
            ax.plot(ddd['q'], ddd['freq'], '-', color=colors[0])
            xlim = [0, 1/3]
            ax.set(
                xticklabels=['$\\Gamma$', 'K'], xticks=xlim, xlim=xlim,
                title=f'AIREBO, AB Stacking'
            )
        ax.plot(r['q'], r['freq'], 'x', color=colors[1])
        ax.set(ylabel='$\\omega~(\\mathrm{cm}^{-1})$')
        fig.tight_layout()
        output = 'band_ab.pdf'
        plt.savefig(f'{output}', bbox_inches='tight')
        if 'png' in output:
            os.system(f'rsub {output}')


# read_pp()
# plot_heat()
# plot_dos()
# plot_band_mo()
plot_band_ab()
