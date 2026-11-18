#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count
import concurrent.futures
from typing import Dict, List
from BLRun.runner import Runner
import os
import pandas as pd

import BLRun as br
yaml.warnings({'YAMLLoadWarning': False})


def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml',
        help='Path to config file')

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config


    with open(config_file, 'r') as conf:
        evaluation = br.ConfigParser.parse(conf)
    print(evaluation)
    print('Evaluation started')


    for idx in range(len(evaluation.runners)):
        evaluation.runners[idx].generateInputs()

    # Sequential execution to avoid distributed computing conflicts
    print(f"Starting {len(evaluation.runners)} algorithms sequentially...")
    
    for idx in range(len(evaluation.runners)):
        algorithm_name = evaluation.runners[idx].name
        print(f"\nAlgorithm {idx+1}/{len(evaluation.runners)}: {algorithm_name}")
        
        proc = evaluation.runners[idx].run()
        if proc:
            print(f"Waiting for {algorithm_name} to complete...")
            proc.wait()
            print(f"{algorithm_name} execution completed.")
        
        print(f"Parsing output for {algorithm_name}...")
        evaluation.runners[idx].parseOutput()
        print(f"{algorithm_name} fully completed.")

    print("\nAll algorithms completed sequentially.")
    print('Evaluation complete')


if __name__ == '__main__':
  main()
