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
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency

# Import optimized utility functions
from .utils import get_truth_edges, get_pred_edges, profile_function

@profile_function
def Spearman(evalObject, algorithmName):
    """
    A function to compute median pairwise Spearman correlation
    of predicted ranked edges, i.e., the outputs of different datasets
    generated from the same reference network, for a given algorithm.
    
    OPTIMIZED VERSION: Uses vectorized pandas operations instead of 
    inefficient permutation-based edge generation and dictionary operations.
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: BLEval
      
    :param algorithmName: Name of the algorithm for which the Spearman correlation is computed.
    :type algorithmName: str
      
    :returns:
        - median: Median of Spearman correlation values
        - mad: Median Absolute Deviation of the Spearman correlation values
    """
    
    print(f"🎯 Computing Spearman correlation for {algorithmName}")
    
    edge_weight_dfs = {}
    
    for dataset in tqdm(evalObject.input_settings.datasets, desc="Processing datasets"):
        dataset_name = dataset['name']
        
        # Handle both single string and list of ground truth networks
        true_edges_files = dataset['trueEdges']
        if isinstance(true_edges_files, str):
            true_edges_files = [true_edges_files]
        elif not isinstance(true_edges_files, list):
            print(f"⚠️  Invalid trueEdges format for {dataset_name}: {type(true_edges_files)}")
            continue
        
        # For Spearman, we use the first/primary ground truth network for gene universe
        # TODO: Consider aggregating across multiple networks in future
        primary_network = true_edges_files[0]
        
        # Load truth edges to get gene universe
        truth_path = str(evalObject.input_settings.datadir) + '/' + dataset_name + '/' + primary_network
        truth_set, truth_df = get_truth_edges(truth_path, directed=True)
        
        print(f"📊 {dataset_name}: {len(truth_set)} true edges")
        
        # Set up algorithm output path
        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset_name + "/" + algorithmName

        rank_path = outDir + "/rankedEdges.csv"
        
        if not os.path.isdir(outDir):
            print(f"⚠️  Output directory {outDir} does not exist for {algorithmName}")
            continue
        
        # Load predicted edges using optimized function
        try:
            pred_df = get_pred_edges(rank_path, directed=True)
            if len(pred_df) == 0:
                print(f"⚠️  No predictions found for {algorithmName} on {dataset_name}")
                continue
        except Exception as e:
            print(f"⚠️  Error loading predictions for {algorithmName} on {dataset_name}: {e}")
            continue
        
        print(f"📈 {dataset_name}: {len(pred_df)} predicted edges loaded")
        
        # Create edge identifier and set as index for fast lookup
        pred_df['EdgeID'] = pred_df['Gene1'] + '|' + pred_df['Gene2']
        pred_indexed = pred_df.set_index('EdgeID')['EdgeWeight']
        
        edge_weight_dfs[dataset_name] = pred_indexed
        print(f"🎯 {dataset_name}: Edge weights indexed for correlation")
    
    if len(edge_weight_dfs) < 2:
        print("⚠️  Need at least 2 datasets for Spearman correlation computation") 
        print("   Single dataset detected - Spearman correlation not applicable")
        return 0.0, 0.0
    
    # Create combined DataFrame with all edge weights
    print("🔄 Combining edge weight data across datasets...")
    combined_df = pd.DataFrame(edge_weight_dfs)
    
    # Fill missing values with 0 (edges not predicted by algorithm get weight 0)
    combined_df = combined_df.fillna(0)
    
    print(f"📊 Combined data: {len(combined_df)} unique edges across {len(combined_df.columns)} datasets")
    
    # Compute pairwise Spearman correlations
    print("🔄 Computing pairwise Spearman correlations...")
    spearman_df = combined_df.corr(method='spearman')
    
    # Extract upper triangular values for statistics
    mask = np.triu(np.ones(spearman_df.shape), k=1).astype(bool)
    correlation_values = spearman_df.values[mask]
    
    # Remove NaN values if any
    correlation_values = correlation_values[~np.isnan(correlation_values)]
    
    if len(correlation_values) == 0:
        print("⚠️  No valid correlations computed")
        return 0.0, 0.0
    
    median_corr = np.median(correlation_values)
    mad_corr = np.median(np.abs(correlation_values - median_corr))
    
    print(f"✅ Spearman correlation - Median: {median_corr:.3f}, MAD: {mad_corr:.3f}")
    print(f"📈 Correlation statistics: min={np.min(correlation_values):.3f}, max={np.max(correlation_values):.3f}")
    
    return median_corr, mad_corr
