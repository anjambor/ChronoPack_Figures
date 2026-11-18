import os
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

# Import optimized utility functions
from .utils import get_truth_edges, get_pred_edges, compute_early_precision_fast, profile_function, resolve_output_path

@profile_function
def EarlyPrec(evalObject, algorithmName, TFEdges = True):
    '''    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the early precision is computed.
    :type algorithmName: str
    
    :param TFEdges: Whether to consider only TF->gene edges
    :type TFEdges: bool
      
    :returns:
        A dictionary containing early precision values
        for a given algorithm for each network configuration.
    '''
    
    print(f"🎯 Computing early precision for {algorithmName}")
    
    EPrecDict = {}
    
    # Use expanded network configurations
    for network_config in tqdm(evalObject.input_settings.network_configs, desc="Processing network configurations"):
        dataset_name = network_config['dataset_name']
        network_file = network_config['network_file']
        column_name = network_config['column_name']
        
        # Load truth edges using optimized function
        truth_path = str(evalObject.input_settings.datadir) + '/' + dataset_name + '/' + network_file
        truth_set, truth_df = get_truth_edges(truth_path, directed=True)
        
        print(f"📊 {column_name}: {len(truth_set)} true edges")
        
        # Set up algorithm output path using smart resolution
        outDir = resolve_output_path(
            str(evalObject.output_settings.base_dir),
            str(evalObject.input_settings.datadir), 
            dataset_name, 
            algorithmName
        )

        rank_path = outDir + "/rankedEdges.csv"
        
        if not os.path.isdir(outDir):
            print(f"⚠️  Output directory {outDir} does not exist for {algorithmName}")
            EPrecDict[column_name] = 0.0
            continue
            
        # Load predicted edges using optimized function
        try:
            pred_df = get_pred_edges(rank_path, directed=True)
            if len(pred_df) == 0:
                print(f"⚠️  No predictions found for {algorithmName} on {column_name}")
                EPrecDict[column_name] = 0.0
                continue
        except Exception as e:
            print(f"⚠️  Error loading predictions for {algorithmName} on {column_name}: {e}")
            EPrecDict[column_name] = 0.0
            continue
        
        print(f"📈 {column_name}: {len(pred_df)} predicted edges loaded")
        
        # Handle TF edges filtering if requested
        if TFEdges:
            tf_genes = set(truth_df['Gene1'].unique())
            print(f"🧬 {column_name}: {len(tf_genes)} TF genes identified")
            
            pred_df = pred_df[pred_df['Gene1'].isin(tf_genes)]
            
            truth_set_tf = {(g1, g2) for g1, g2 in truth_set if g1 in tf_genes}
            truth_set = truth_set_tf
            
            print(f"🎯 {column_name}: {len(truth_set)} TF->gene true edges, {len(pred_df)} TF->gene predictions")
        
        k = len(truth_set)
        result = compute_early_precision_fast(truth_set, pred_df, k=k)
        EPrecDict[column_name] = result['precision']

        print(f"✅ {column_name}: Early Precision = {result['precision']:.3f}, Early Recall = {result['recall']:.3f}")
    
    return EPrecDict
