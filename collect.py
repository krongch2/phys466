import os
import re

def read_fc():

    workdir = os.getcwd()
    os.chdir('data')
    for dirname in sorted(os.listdir('.')):

        os.chdir(dirname)
        if os.path.isfile('FORCE_CONSTANTS') and not os.path.isfile('pp.pkl'):
            print(dirname)
            with open('pp.slurm', 'w') as f:
                f.write(
f'''#! /bin/bash
#SBATCH --job-name="{dirname}"
#SBATCH --time=4:00:00
#SBATCH --partition="secondary"
#SBATCH --cpus-per-task=20

cd $SLURM_SUBMIT_DIR

python {workdir}/process.py
''')
            # os.system('sbatch pp.slurm')
            os.system(f'python {workdir}/process.py')

        os.chdir('..')
    os.chdir('..')

read_fc()
