#!/usr/bin/env python3
'''
Here we provide a script to launch Vacuum Emission simulation
IN PARALLEL (using submitit for semi-slurm submission type),
do postprocessing and measure performance
'''
import argparse

import submitit

from quvac.simulation import quvac_simulation


def parse_args():
    description = "Calculate quantum vacuum signal for given external fields"
    argparser = argparse.ArgumentParser(description=description)
    argparser.add_argument("--input", "-i", default=None,
                           help="Input yaml file with field and grid params")
    argparser.add_argument("--output", "-o", default=None,
                           help="Path to save simulation data to")
    argparser.add_argument("--wisdom", default='wisdom/fftw-wisdom',
                           help="File to save pyfftw-wisdom")
    return argparser.parse_args()


def quvac_sim_manager(ini_file, save_path=None, wisdom_file='wisdom/fftw-wisdom'):
    '''
    We read ini script and parallelization parameters
    Depending on available jobs, we split the total time interval
    into several sub-intervals and submit each sub-interval for 
    calculation as a separate quvac simulation (without postprocessing).
    Then we gather all S1 and S2 in the main process and do main postprocessing.
    '''
    


if __name__ == '__main__':
    args = parse_args()
    quvac_simulation(args.input, args.output, args.wisdom)


