import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rc('font', family='serif')
plt.rc('text', usetex=True)
colors = sns.color_palette()

thz_to_cm1 = 33.35641

def plot_therm(csv, output):
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

# if os.path.isfile('tp.pkl'):
        #     with open('tp.pkl', 'rb') as f:
        #         tp_dict = pickle.load(f)

        #     with open('POSCAR_unitcell') as f:
        #         text = f.read()
        #         for line in text.split('\n'):
        #             m =  re.search('^\d+$', line)
        #             if m:
        #                 natoms = m.group(0)
        #     pot, angle = dirname.split('_')
        #     d = pd.DataFrame(tp_dict)
        #     d['pot'] = pot
        #     d['angle'] = float(angle)
        #     d['natoms'] = int(natoms)
        #     l.append(d)


# def myband():
    #     d = phonon.get_band_structure_dict()
    #     q = d['qpoints']
    #     freq = d['frequencies']
    #     fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))
    #     x = q[-1]
    #     y = freq[-1]
    #     for i in range(y.shape[1]):
    #         ax.plot(x, y[:, i]*thz_to_cm1, color=colors[0])
    #     xlim = [0, 1/3]
    #     ax.set(
    #         xticklabels=labels, xticks=xlim, xlim=xlim,
    #         ylabel='Phonon frequency~$\\mathrm{cm}^{-1}$'
    #     )
    #     fig.tight_layout()
    #     plt.savefig('myband.png', bbox_inches='tight')
    #     os.system('rsub myband.png')