#!/usr/bin/env python
# coding: utf-8

import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

# local imports
import BLEval as ev 
from BLEval.utils import profile_function

def get_parser() -> argparse.ArgumentParser:
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('-c','--config', default='config.yaml',
        help="Configuration file containing list of datasets "
              "algorithms and output specifications.\n")
    
    parser.add_argument('-a', '--auc', action="store_true", default=False,
        help="Compute median of areas under Precision-Recall and ROC curves.\n")
    
    parser.add_argument('-j', '--jaccard', action="store_true", default=False,
      help="Compute median Jaccard index of predicted top-k networks "
      "for each algorithm for a given set of datasets generated "
      "from the same ground truth network.\n")

    parser.add_argument('-r', '--spearman', action="store_true", default=False,
      help="Compute median Spearman Corr. of predicted edges "
      "for each algorithm  for a given set of datasets generated "
      " from the same ground truth network.\n")

    parser.add_argument('-t', '--time', action="store_true", default=False,
      help="Analyze time taken by each algorithm for a.\n")
    
    parser.add_argument('-e', '--epr', action="store_true", default=False,
      help="Compute median early precision.")
    
    parser.add_argument('-s','--sepr', action="store_true", default=False,
      help="Analyze median (signed) early precision for activation and inhibitory edges.")

    parser.add_argument('-m','--motifs', action="store_true", default=False,
      help="Compute network motifs in the predicted top-k networks.")
    
    parser.add_argument('-p','--paths', action="store_true", default=False,
      help="Compute path length statistics on the predicted top-k networks.")

    parser.add_argument('-b','--borda', action="store_true", default=False,
      help="Compute edge ranked list using the various Borda aggregatio methods.")

    parser.add_argument('--parallel', action="store_true", default=False,
      help="Enable parallel processing for evaluation metrics.")
    
    parser.add_argument('--max-workers', type=int, default=None,
      help="Maximum number of parallel workers (default: number of CPU cores).")
        
    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()
    
    return opts

@profile_function 
def compute_metric_parallel(metric_func, evalSummarizer, *args, **kwargs):
    """Wrapper function for parallel metric computation."""
    return metric_func(*args, **kwargs)

@profile_function
def main():
    opts = parse_arguments()
    config_file = opts.config

    evalConfig = None

    with open(config_file, 'r') as conf:
        evalConfig = ev.ConfigParser.parse(conf)
        
    print('\n🚀 Post-run evaluation started...')
    print(f'📊 Evaluating {len(evalConfig.input_settings.datasets)} datasets with {len(evalConfig.input_settings.algorithms)} algorithms')
    
    evalSummarizer = ev.BLEval(evalConfig.input_settings, evalConfig.output_settings)
    
    outDir = str(evalSummarizer.output_settings.base_dir) + \
            str(evalSummarizer.input_settings.datadir).split("inputs")[1] + "/"+\
            str(evalSummarizer.output_settings.output_prefix) + "-"
    
    # Ensure output directory exists
    output_dir_path = Path(outDir).parent
    output_dir_path.mkdir(parents=True, exist_ok=True)
    print(f"📁 Output directory: {output_dir_path}")
    print(f"📝 Output file prefix: {outDir}")
    
    # Collect all requested metrics for potential parallel execution
    metrics_to_compute = []
    
    if opts.auc:
        metrics_to_compute.append(('AUC', evalSummarizer.computeAUC, 'AUPRC', 'AUROC'))

    if opts.epr:
        metrics_to_compute.append(('EarlyPrecision', evalSummarizer.computeEarlyPrec, 'EPr', None))
    
    if opts.jaccard:
        metrics_to_compute.append(('Jaccard', evalSummarizer.computeJaccard, 'Jaccard', None))
        
    if opts.spearman:
        metrics_to_compute.append(('Spearman', evalSummarizer.computeSpearman, 'Spearman', None))
        

    # Execute metrics (parallel or sequential)
    if opts.parallel and len(metrics_to_compute) > 1:
        print(f"🔄 Running {len(metrics_to_compute)} metrics in parallel...")
        
        max_workers = opts.max_workers or min(len(metrics_to_compute), cpu_count())
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_metric = {}
            
            for metric_name, metric_func, output1, output2 in metrics_to_compute:
                future = executor.submit(compute_metric_parallel, metric_func, evalSummarizer)
                future_to_metric[future] = (metric_name, output1, output2)
            
            # Collect results
            for future in concurrent.futures.as_completed(future_to_metric):
                metric_name, output1, output2 = future_to_metric[future]
                try:
                    result = future.result()
                    
                    # Save results
                    if output2:  # Two outputs (like AUC)
                        result[0].to_csv(outDir + output1 + '.csv')
                        result[1].to_csv(outDir + output2 + '.csv')
                        print(f"✅ {metric_name} completed - saved {output1}.csv and {output2}.csv")
                    else:  # Single output
                        result.to_csv(outDir + output1 + '.csv')
                        print(f"✅ {metric_name} completed - saved {output1}.csv")
                        
                except Exception as exc:
                    print(f"❌ {metric_name} generated an exception: {exc}")
    else:
        # Sequential execution (original behavior)
        print("🔄 Running metrics sequentially...")
        
        # Compute and plot ROC, PRC and report median AUROC, AUPRC    
        if opts.auc:
            print('\n\n📈 Computing areas under ROC and PR curves...')
            AUPRC, AUROC = evalSummarizer.computeAUC()
            AUPRC.to_csv(outDir+'AUPRC.csv')
            AUROC.to_csv(outDir+'AUROC.csv')
            print('✅ AUC computation completed')
        
        # Compute early precision
        if opts.epr:
            print('\n\n🎯 Computing early precision values...')
            ePRDF = evalSummarizer.computeEarlyPrec()
            ePRDF.to_csv(outDir + "EPr.csv")
            print('✅ Early precision computation completed')
            
        # Compute Jaccard index    
        if opts.jaccard:
            print('\n\n🎯 Computing Jaccard index...')
            jaccDict = evalSummarizer.computeJaccard()
            jaccDict.to_csv(outDir + "Jaccard.csv")
            print('✅ Jaccard computation completed')
            
        # Compute Spearman correlation scores
        if opts.spearman:
            print('\n\n📊 Computing Spearman\'s correlation...')
            corrDict = evalSummarizer.computeSpearman()
            corrDict.to_csv(outDir + "Spearman.csv")
            print('✅ Spearman computation completed')
        
    
    # Handle remaining metrics that don't benefit from parallelism as much
    
    # Compute median time taken
    if opts.time:
        print('\n\n⏱️  Computing time taken...')
        TimeDict = evalSummarizer.parseTime()
        pd.DataFrame(TimeDict).to_csv(outDir+'Times.csv')
        print('✅ Time analysis completed')
                        
    # Compute early precision for activation and inhibitory edges
    if opts.sepr:
        print('\n\n🧬 Computing early precision values for activation and inhibitory edges...')
        actDF, inhDF = evalSummarizer.computeSignedEPrec()
        actDF.to_csv(outDir + "EPr-Activation.csv")
        inhDF.to_csv(outDir + "EPr-Inhibitory.csv")
        print('✅ Signed early precision computation completed')

    # Compute median time taken
    if opts.motifs:
        print('\n\n🔗 Computing network motifs...')
        FBL, FFL, MI = evalSummarizer.computeNetMotifs()
        FBL.to_csv(outDir+'NetworkMotifs-FBL.csv')
        FFL.to_csv(outDir+'NetworkMotifs-FFL.csv')
        MI.to_csv(outDir+'NetworkMotifs-MI.csv')
        print('✅ Network motifs computation completed')
        
    # Compute path statistics such as number of TP, FP, 
    # and path lengths among TP in the top-k networks.
    if opts.paths:
        print('\n\n🛤️  Computing path length statistics on predicted networks...')
        evalSummarizer.computePaths()
        print('✅ Path statistics computation completed')
        
   # Compute edge ranked list using the borda method
    if opts.borda:
        print('\n\n🗳️  Computing edge ranked list using the borda method')
        evalSummarizer.computeBorda()
        print('✅ Borda aggregation completed')

    print('\n\n🎉 Evaluation complete! Results saved to:', outDir + '*')
    print('📁 Check the output directory for CSV files with evaluation metrics.')


if __name__ == '__main__':
  main()
