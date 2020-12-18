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

{min_section}
'''

pots = {
'lj':
'''pair_style lj/cut 4.0
pair_coeff * * 1.0 1.0 4.0''',

'rebo':
'''pair_style hybrid/overlay rebo kolmogorov/crespi/full 16.0 1
pair_coeff * * rebo {pot_path}/CH.rebo C # chemical
pair_coeff * * kolmogorov/crespi/full {pot_path}/CH_taper.KC C # long-range''',

'airebo':
'''pair_style airebo 2.0
pair_coeff * * {pot_path}/CH.airebo C''',

'tersoff':
'''pair_style   tersoff
pair_coeff   * * {pot_path}/SiC.tersoff C'''
}

class Job(object):

    def __init__(self, options):
        self.geom_path = '/home/krongch2/projects/phys466/geoms'
        self.pot_path = '/home/krongch2/projects/phys466/potentials'
        self.data_path = 'data'
        self.angle = 0
        self.min_section = ''

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
        self.geom_tail = f'{self.angle:04.1f}.pos'
        self.geom_full_path = os.path.join(self.geom_path, self.geom_tail)
        self.dirname = os.path.join(self.data_path, f'{self.pot}_{self.angle:04.1f}')
        self.pot_section = pots[self.pot].format(**self.__dict__)

    def run_min(self):
        workdir = os.getcwd()
        os.makedirs(self.dirname, exist_ok=True)
        os.chdir(self.dirname)
        self.geom = 'input.pos'
        os.system(f'cp {self.geom_full_path} {self.geom}')
        self.min_section = '''reset_timestep 0
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10
thermo_style multi

min_style sd
minimize 1.0e-4 1.0e-6 1000 10000

write_data relax.pos nocoeff
'''

        self.min_section = '''fix 1 all box/relax iso 0.0 vmax 0.001

# Timestep to make the integration of the motion equation
timestep        0.0005

# Parameters to print out on the screen and log.lammps file
thermo_style    custom step temp etotal vol lx ly lz press pxx pyy pzz cpu
thermo          1000

group top type 1
group bot type 2

# Energy minimization parameters
min_style       cg
minimize        1e-15 1e-12 100000 10000

unfix 1
min_style       fire
minimize        1e-15 1e-12 100000 10000

fix 1 all box/relax iso 0.0 vmax 0.001
min_style       cg
minimize        1e-15 1e-12 100000 10000

unfix 1
min_style       fire
minimize        1e-15 1e-12 100000 10000

write_data relax.pos nocoeff

variable Lx equal "lx"
variable Ly equal "ly"

variable teng equal "pe"
variable natoms equal "count(all)"
variable ecoh equal "v_teng/v_natoms"
# print "${teng}, ${ecoh}" file energy.dat
# print "${Lx}, ${Ly}" file Lxy.dat

        '''
        with open('relax.in', 'w') as f:
            f.write(template.format(**self.__dict__))

        os.system('lmp_mpi -in relax.in')
        os.chdir(workdir)

    def run_phlmp(self):
        '''
        Runs phonoLAMMPS
        '''
        workdir = os.getcwd()
        os.makedirs(self.dirname, exist_ok=True)
        os.chdir(self.dirname)

        self.min_section = ''
        self.geom = 'relax.pos'

        with open('phlmp.in', 'w') as f:
            f.write(template.format(**self.__dict__))

        os.system('lmp_mpi -in phlmp.in')

        with open('phlmp.slurm', 'w') as f:
            f.write(
f'''#! /bin/bash
#SBATCH --job-name="{self.pot}_{self.angle}"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20
#SBATCH -e %x.e%A

cd $SLURM_SUBMIT_DIR

phonolammps phlmp.in --dim 2 2 2 -c POSCAR_unitcell
''')
        os.system('phonolammps phlmp.in --dim 2 2 2 -c POSCAR_unitcell')
        # os.system('sbatch phlmp.slurm')
        os.chdir(workdir)

params = {
    'pot': ['lj', 'airebo', 'tersoff'],
    # 'pot': ['lj', 'airebo'],
    # 'angle': [0.0, 3.5, 4.4, 7.3, 21.7, 13.2]
    'angle': [0.0, 13.2, 21.7]
}

for p in itertools.product(*params.values()):
    options = dict(zip(params.keys(), p))
    job = Job(options)

    if os.path.isfile(os.path.join(job.dirname, 'FORCE_CONSTANTS')):
        continue
    job.run_min()
    job.run_phlmp()
