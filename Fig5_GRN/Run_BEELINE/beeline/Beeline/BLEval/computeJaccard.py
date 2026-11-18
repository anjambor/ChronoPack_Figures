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
from .utils import get_truth_edges, get_pred_edges, compute_jaccard_fast, profile_function

@profile_function
def Jaccard(evalObject, algorithmName):
    """
    A function to compute median pairwise Jaccard similarity index
    of predicted top-k edges for a given set of datasets (obtained from
    the same reference network). Here k is the number of edges in the
    reference network (excluding self loops). 
    
    OPTIMIZED VERSION: Uses fast set operations instead of 
    inefficient permutation-based edge generation.
    
    :param evalObject: An object of class :class:`BLEval.BLEval`.
    :type evalObject: :obj:`BLEval`
      
    :param algorithmName: Name of the algorithm for which the Jaccard similarity is computed.
    :type algorithmName: str
      
    :returns:
        - median: Median of Jaccard correlation values
        - mad: Median Absolute Deviation of the Jaccard correlation values
    """
    
    print(f"🎯 Computing Jaccard similarity for {algorithmName}")
    
    pred_sets = {}
    
    for dataset in tqdm(evalObject.input_settings.datasets, desc="Processing datasets"):
        dataset_name = dataset['name']
        
        # Handle both single string and list of ground truth networks
        true_edges_files = dataset['trueEdges']
        if isinstance(true_edges_files, str):
            true_edges_files = [true_edges_files]
        elif not isinstance(true_edges_files, list):
            print(f"⚠️  Invalid trueEdges format for {dataset_name}: {type(true_edges_files)}")
            continue
        
        # For Jaccard, we typically use the first/primary ground truth network
        # TODO: Consider aggregating across multiple networks in future
        primary_network = true_edges_files[0]
        
        # Load truth edges using optimized function
        truth_path = str(evalObject.input_settings.datadir) + '/' + dataset_name + '/' + primary_network
        truth_set, truth_df = get_truth_edges(truth_path, directed=True)
        
        num_true_edges = len(truth_set)
        print(f"📊 {dataset_name}: {num_true_edges} true edges")
        
        # Set up algorithm output path
        outDir = str(evalObject.output_settings.base_dir) + \
                 str(evalObject.input_settings.datadir).split("inputs")[1] + \
                 "/" + dataset_name + "/" + algorithmName

        rank_path = outDir + "/rankedEdges.csv"
        
        if not os.path.isdir(outDir):
            print(f"⚠️  Output directory {outDir} does not exist for {algorithmName}")
            pred_sets[dataset_name] = set()
            continue
        
        # Load predicted edges using optimized function
        try:
            pred_df = get_pred_edges(rank_path, directed=True)
            if len(pred_df) == 0:
                print(f"⚠️  No predictions found for {algorithmName} on {dataset_name}")
                pred_sets[dataset_name] = set()
                continue
        except Exception as e:
            print(f"⚠️  Error loading predictions for {algorithmName} on {dataset_name}: {e}")
            pred_sets[dataset_name] = set()
            continue
        
        print(f"📈 {dataset_name}: {len(pred_df)} predicted edges loaded")
        
        # Get top-k predictions where k = number of true edges
        k = min(len(pred_df), num_true_edges)
        
        if k > 0:
            # Get weight threshold for top-k edges
            weight_threshold = pred_df.iloc[k-1]['EdgeWeight']
            
            # Include all edges with weight >= threshold (handles ties)
            top_k_df = pred_df.loc[pred_df['EdgeWeight'] >= weight_threshold]
            
            # Convert to set of edge tuples
            pred_set = set(zip(top_k_df['Gene1'], top_k_df['Gene2']))
            pred_sets[dataset_name] = pred_set
            
            print(f"🎯 {dataset_name}: Top-{k} edges selected ({len(pred_set)} edges including ties)")
        else:
            pred_sets[dataset_name] = set()
    
    # Compute pairwise Jaccard similarities using optimized function
    if len(pred_sets) < 2:
        print("⚠️  Need at least 2 datasets for Jaccard similarity computation")
        print("   Single dataset detected - Jaccard similarity not applicable")
        return 0.0, 0.0
    
    jaccard_df = compute_jaccard_fast(pred_sets)
    
    # Extract upper triangular values for statistics
    mask = np.triu(np.ones(jaccard_df.shape), k=1).astype(bool)
    jaccard_values = jaccard_df.values[mask]
    
    if len(jaccard_values) == 0:
        return 0.0, 0.0
    
    median_jaccard = np.median(jaccard_values)
    mad_jaccard = np.median(np.abs(jaccard_values - median_jaccard))
    
    print(f"✅ Jaccard similarity - Median: {median_jaccard:.3f}, MAD: {mad_jaccard:.3f}")
    
    return median_jaccard, mad_jaccard

def computePairwiseJacc(inDict):
    """
    DEPRECATED: Use compute_jaccard_fast from utils instead.
    
    A helper function to compute all pairwise Jaccard similarity indices
    of predicted top-k edges for a given set of datasets (obtained from
    the same reference network). Here k is the number of edges in the
    reference network (excluding self loops). 
    
    :param inDict:  A dictionary containing top-k predicted edges for each dataset. Here, keys are the dataset name and the values are the set of top-k edges.
    :type inDict: dict
    :returns:
        A dataframe containing pairwise Jaccard similarity index values
    """
    print("⚠️  WARNING: Using deprecated computePairwiseJacc function.")
    print("   Use compute_jaccard_fast from utils instead for better performance.")
    
    # Redirect to optimized version
    from .utils import compute_jaccard_fast
    return compute_jaccard_fast(inDict)
