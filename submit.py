import itertools
import os


template = '''
# 3D copper block simulation
boundary     p p p
units        metal
atom_style   atomic

# geometry
read_data    {geom}

mass        *   12.0107

# EAM potential
{pot_section}
neighbor     2.0 nsq
neigh_modify every 2 delay 10 check yes
'''

pots = {
'lj':
'''pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
''',

'rebo':
'''pair_style hybrid/overlay rebo kolmogorov/crespi/full 16.0 1
pair_coeff * * rebo {pot_path}/CH.rebo C # chemical
pair_coeff * * kolmogorov/crespi/full {pot_path}/CH_taper.KC C # long-range
''',

'tersoff':
'''pair_style   tersoff
pair_coeff   * * {pot_path}/SiC.tersoff C
'''
}

class Job(object):

    def __init__(self, options):
        self.geom_path = '/home/krongch2/projects/phys466/geoms'
        self.pot_path = '/home/krongch2/projects/phys466/potentials'
        self.data_path = 'data'
        self.angle = 0

        self.set_options(options)
        self.predefine()

    def set_options(self, options):
        '''
        Sets the attributes from options dict
        '''
        def to_numeric(x):
            try:
                v = float(x)
                if '.' not in x:
                    v = int(x)
            except:
                v = x
            return v

        for k in options.keys():
            if not k in self.__dict__.keys():
                pass
                # print(f'"{k}" is not a keyword for Job')
                # raise AssertionError(f'"{k}" is not a keyword for Job')

            setattr(self, k, to_numeric(options[k]))

    def predefine(self):
        '''
        Sets dependent attributes from the starting options
        '''
        self.geom = os.path.join(self.geom_path, f'{self.angle:04.1f}.pos')
        self.dirname = os.path.join(self.data_path, f'{self.pot}_{self.angle:04.1f}')
        self.pot_section = pots[self.pot].format(**self.__dict__)

    def run_phlmp(self):
        '''
        Runs phonoLAMMPS
        '''
        workdir = os.getcwd()
        os.makedirs(self.dirname, exist_ok=True)
        os.chdir(self.dirname)

        with open('lmp.in', 'w') as f:
            f.write(template.format(**self.__dict__))

        os.system('lmp_mpi -in lmp.in')

        os.chdir(workdir)

# def run_phlmp():

#     with open('phlmp.slurm', 'w') as f:
#         f.write(
# '''#! /bin/bash
# #SBATCH --job-name="phlmp"
# #SBATCH --time=4:00:00
# #SBATCH --partition="secondary"
# #SBATCH --cpus-per-task=1

# cd $SLURM_SUBMIT_DIR

# ''')

params = {
    'pot': ['lj', 'rebo', 'tersoff'],
    'angle': [3.5, 4.4, 7.3, 21.7, 31.2]
}

for p in itertools.product(*params.values()):
    options = dict(zip(params.keys(), p))
    job = Job(options)
    job.run_phlmp()
