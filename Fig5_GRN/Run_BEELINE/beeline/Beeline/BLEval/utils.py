import pandas as pd
import numpy as np
from pathlib import Path
from typing import Set, Tuple, Dict, Optional, Union, List
import time
import cProfile
import pstats
from functools import wraps

def profile_function(func):
    """Decorator to profile function execution time and memory usage."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        
        # Run with profiling
        profiler = cProfile.Profile()
        profiler.enable()
        result = func(*args, **kwargs)
        profiler.disable()
        
        end_time = time.time()
        
        # Print timing information
        print(f"\n⏱️  {func.__name__} completed in {end_time - start_time:.2f} seconds")
        
        # Optional: save detailed profiling stats
        stats = pstats.Stats(profiler)
        stats.sort_stats('cumulative')
        # Uncomment to see detailed profiling:
        # stats.print_stats(10)
        
        return result
    return wrapper

def resolve_output_path(base_dir: str, datadir: str, dataset_name: str, algorithm_name: str) -> str:
    """
    Smart path resolution for algorithm output directories.
    
    :param base_dir: Base output directory (e.g., 'outputs')
    :param datadir: Input data directory (e.g., 'inputs/scRNA-Seq')
    :param dataset_name: Dataset name (e.g., 'K562')
    :param algorithm_name: Algorithm name (e.g., 'GENIE3')
    :returns: Path to algorithm output directory
    """
    # Extract the suffix after 'inputs'
    if "inputs/" in datadir:
        datadir_suffix = datadir.split("inputs/")[1]
    elif "inputs" in datadir:
        datadir_suffix = datadir.split("inputs")[1].lstrip("/")
    else:
        datadir_suffix = datadir
    
    # Try different path combinations
    path_options = [
        f"{base_dir}/{datadir_suffix}/{dataset_name}/{algorithm_name}",  # Direct path
        f"../{base_dir}/{datadir_suffix}/{dataset_name}/{algorithm_name}",  # Parent directory
        f"outputs/{datadir_suffix}/{dataset_name}/{algorithm_name}",  # Hardcoded outputs
        f"../outputs/{datadir_suffix}/{dataset_name}/{algorithm_name}"  # Parent + hardcoded
    ]
    
    for path_option in path_options:
        if Path(path_option).exists():
            return path_option
    
    # If none exist, return the first option (most likely correct)
    return path_options[0]

def get_truth_edges(truth_path: Union[str, Path], directed: bool = True) -> Tuple[Set[Tuple[str, str]], pd.DataFrame]:
    """
    Load truth edges from CSV file and return both set and DataFrame representations.
    
    :param truth_path: Path to the truth edges CSV file
    :param directed: Whether to treat edges as directed
    :returns: Tuple of (edge_set, truth_dataframe)
    """
    truth_df = pd.read_csv(truth_path, sep=',', header=0, index_col=None)
    
    # Convert gene names to uppercase for case-insensitive matching
    truth_df['Gene1'] = truth_df['Gene1'].astype(str).str.upper()
    truth_df['Gene2'] = truth_df['Gene2'].astype(str).str.upper()
    
    # Remove self-edges and duplicates
    truth_df = truth_df.loc[truth_df['Gene1'] != truth_df['Gene2']]
    truth_df.drop_duplicates(keep='first', inplace=True)
    truth_df.reset_index(drop=True, inplace=True)
    
    if directed:
        # For directed graphs, keep edges as (Gene1, Gene2)
        edge_set = set(zip(truth_df['Gene1'], truth_df['Gene2']))
    else:
        # For undirected graphs, normalize edge representation
        edge_set = set()
        for _, row in truth_df.iterrows():
            edge = tuple(sorted([row['Gene1'], row['Gene2']]))
            edge_set.add(edge)
    
    # Create DataFrame with y_true labels
    if directed:
        truth_labeled = truth_df.copy()
    else:
        # For undirected, create normalized representation
        truth_labeled = pd.DataFrame([
            {'Gene1': min(g1, g2), 'Gene2': max(g1, g2)} 
            for g1, g2 in edge_set
        ])
    
    truth_labeled['y_true'] = 1
    
    return edge_set, truth_labeled

def get_pred_edges(pred_path: Union[str, Path], top_k: Optional[int] = None, 
                   directed: bool = True) -> pd.DataFrame:
    """
    Load predicted edges from CSV file with optional top-k filtering.
    
    :param pred_path: Path to the predicted edges CSV file
    :param top_k: Maximum number of top edges to keep (None for all)
    :param directed: Whether to treat edges as directed
    :returns: DataFrame with Gene1, Gene2, EdgeWeight columns
    """
    if not Path(pred_path).exists():
        return pd.DataFrame(columns=['Gene1', 'Gene2', 'EdgeWeight'])
    
    pred_df = pd.read_csv(pred_path, sep='\t', header=0, index_col=None)
    
    # Convert gene names to uppercase for case-insensitive matching
    pred_df['Gene1'] = pred_df['Gene1'].astype(str).str.upper()
    pred_df['Gene2'] = pred_df['Gene2'].astype(str).str.upper()
    
    # Remove self-edges and duplicates
    pred_df = pred_df.loc[pred_df['Gene1'] != pred_df['Gene2']]
    pred_df.drop_duplicates(keep='first', inplace=True)
    
    # Clean edge weights
    pred_df['EdgeWeight'] = pred_df['EdgeWeight'].abs().round(6)
    
    # Sort by edge weight (descending)
    pred_df = pred_df.sort_values('EdgeWeight', ascending=False)
    pred_df.reset_index(drop=True, inplace=True)
    
    # Apply top-k filtering if specified
    if top_k is not None and len(pred_df) > top_k:
        # Find the weight threshold for top-k edges
        weight_threshold = pred_df.iloc[top_k-1]['EdgeWeight']
        # Include all edges with weight >= threshold (handles ties)
        pred_df = pred_df.loc[pred_df['EdgeWeight'] >= weight_threshold]
    
    if not directed:
        # For undirected graphs, normalize edge representation
        pred_normalized = []
        for _, row in pred_df.iterrows():
            gene1, gene2 = sorted([row['Gene1'], row['Gene2']])
            pred_normalized.append({
                'Gene1': gene1,
                'Gene2': gene2, 
                'EdgeWeight': row['EdgeWeight']
            })
        pred_df = pd.DataFrame(pred_normalized)
        # Remove duplicates that may have been created by normalization, keeping max weight
        pred_df = pred_df.groupby(['Gene1', 'Gene2'])['EdgeWeight'].max().reset_index()
    
    return pred_df

def make_truth_pred_df(truth_df: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge truth and prediction DataFrames for sklearn evaluation.
    
    :param truth_df: DataFrame with Gene1, Gene2, y_true columns
    :param pred_df: DataFrame with Gene1, Gene2, EdgeWeight columns
    :returns: DataFrame with y_true and y_score columns ready for sklearn
    """
    # Merge truth and predictions
    merged = truth_df.merge(pred_df, on=['Gene1', 'Gene2'], how='outer')
    
    # Fill missing values
    merged['y_true'] = merged['y_true'].fillna(0).astype(int)
    merged['y_score'] = merged['EdgeWeight'].fillna(0.0)
    
    return merged[['y_true', 'y_score']]

def compute_early_precision_fast(truth_set: Set[Tuple[str, str]], 
                                pred_df: pd.DataFrame, 
                                k: Optional[int] = None) -> Dict[str, float]:
    """
    Fast computation of early precision using set operations.
    
    :param truth_set: Set of true edges
    :param pred_df: DataFrame with ranked predictions
    :param k: Number of top predictions to consider (default: len(truth_set))
    :returns: Dictionary with precision and recall values
    """
    if len(pred_df) == 0:
        return {'precision': 0.0, 'recall': 0.0}
    
    if k is None:
        k = len(truth_set)
    
    # Get top-k predictions
    top_k_df = pred_df.head(k)
    pred_set = set(zip(top_k_df['Gene1'], top_k_df['Gene2']))
    
    # Compute intersection
    intersection = truth_set.intersection(pred_set)
    
    precision = len(intersection) / len(pred_set) if len(pred_set) > 0 else 0.0
    recall = len(intersection) / len(truth_set) if len(truth_set) > 0 else 0.0
    
    return {'precision': precision, 'recall': recall}

def compute_jaccard_fast(pred_sets: Dict[str, Set[Tuple[str, str]]]) -> pd.DataFrame:
    """
    Fast computation of pairwise Jaccard similarities using set operations.
    
    :param pred_sets: Dictionary mapping dataset names to edge sets
    :returns: DataFrame with pairwise Jaccard similarities
    """
    dataset_names = list(pred_sets.keys())
    n_datasets = len(dataset_names)
    
    jaccard_matrix = np.zeros((n_datasets, n_datasets))
    
    for i, name_i in enumerate(dataset_names):
        for j, name_j in enumerate(dataset_names):
            set_i, set_j = pred_sets[name_i], pred_sets[name_j]
            
            if len(set_i) == 0 and len(set_j) == 0:
                jaccard_matrix[i, j] = 1.0  # Both empty sets
            elif len(set_i) == 0 or len(set_j) == 0:
                jaccard_matrix[i, j] = 0.0  # One empty, one non-empty
            else:
                intersection = len(set_i.intersection(set_j))
                union = len(set_i.union(set_j))
                jaccard_matrix[i, j] = intersection / union
    
    return pd.DataFrame(jaccard_matrix, index=dataset_names, columns=dataset_names)

def get_unique_genes(truth_path: Union[str, Path]) -> np.ndarray:
    """
    Get unique genes from truth edges file.
    
    :param truth_path: Path to truth edges CSV
    :returns: Array of unique gene names (uppercase)
    """
    truth_df = pd.read_csv(truth_path, sep=',', header=0, index_col=None)
    
    # Convert gene names to uppercase for case-insensitive matching
    truth_df['Gene1'] = truth_df['Gene1'].astype(str).str.upper()
    truth_df['Gene2'] = truth_df['Gene2'].astype(str).str.upper()
    
    return np.unique(truth_df[['Gene1', 'Gene2']].values)

def log_performance_improvement(func_name: str, old_time: float, new_time: float, 
                              old_memory: Optional[float] = None, new_memory: Optional[float] = None):
    """
    Log performance improvements for optimization tracking.
    """
    speedup = old_time / new_time if new_time > 0 else float('inf')
    print(f"\n🚀 Performance Improvement for {func_name}:")
    print(f"   ⏱️  Time: {old_time:.2f}s → {new_time:.2f}s ({speedup:.1f}x speedup)")
    
    if old_memory and new_memory:
        memory_reduction = (old_memory - new_memory) / old_memory * 100
        print(f"   🧠 Memory: {old_memory:.1f}MB → {new_memory:.1f}MB ({memory_reduction:.1f}% reduction)")

# ====================================================================
# NEW FUNCTIONS: Multiple Ground Truth Network Support
# ====================================================================

def extract_network_name(filename: str) -> str:
    """
    Extract network name from filename by removing file extension.
    
    :param filename: Network filename (e.g., "ChIP-seq-network.csv")
    :returns: Network name without extension (e.g., "ChIP-seq-network")
    """
    # Simply return filename without extension
    return Path(filename).stem

def generate_column_name(dataset_name: str, network_name: str) -> str:
    """
    Generate consistent column names for evaluation results.
    
    :param dataset_name: Dataset name (e.g., "K562")
    :param network_name: Network name (e.g., "ChIP-seq")
    :returns: Column name (e.g., "K562_ChIP-seq")
    """
    return f"{dataset_name}_{network_name}"

def parse_multiple_networks(dataset_config: Dict) -> List[Dict]:
    """
    Parse trueEdges as either string or list and return list of network configurations.
    
    :param dataset_config: Dataset configuration from YAML
    :returns: List of network configurations with individual files
    """
    true_edges = dataset_config.get('trueEdges')
    
    if isinstance(true_edges, str):
        # Single network - backward compatibility
        network_name = extract_network_name(true_edges)
        return [{
            'dataset_name': dataset_config['name'],
            'network_name': network_name,
            'network_file': true_edges,
            'column_name': generate_column_name(dataset_config['name'], network_name)
        }]
    
    elif isinstance(true_edges, list):
        # Multiple networks
        networks = []
        for network_file in true_edges:
            network_name = extract_network_name(network_file)
            networks.append({
                'dataset_name': dataset_config['name'],
                'network_name': network_name, 
                'network_file': network_file,
                'column_name': generate_column_name(dataset_config['name'], network_name)
            })
        return networks
    
    else:
        raise ValueError(f"trueEdges must be string or list, got {type(true_edges)}")

def get_multiple_truth_edges(datadir: Union[str, Path], dataset_name: str, 
                           network_files: Union[str, List[str]], 
                           directed: bool = True) -> Dict[str, Tuple[Set[Tuple[str, str]], pd.DataFrame]]:
    """
    Load multiple truth edge networks for a dataset.
    
    :param datadir: Base data directory
    :param dataset_name: Dataset name
    :param network_files: Single file or list of network files
    :param directed: Whether to treat edges as directed
    :returns: Dictionary mapping network names to (edge_set, dataframe) tuples
    """
    if isinstance(network_files, str):
        network_files = [network_files]
    
    results = {}
    
    for network_file in network_files:
        network_name = extract_network_name(network_file)
        truth_path = Path(datadir) / dataset_name / network_file
        
        try:
            edge_set, truth_df = get_truth_edges(truth_path, directed=directed)
            results[network_name] = (edge_set, truth_df)
            print(f"📊 Loaded {len(edge_set)} edges from {network_name} network")
        except Exception as e:
            print(f"⚠️  Failed to load {network_name} network from {truth_path}: {e}")
            # Add empty network to maintain consistency
            results[network_name] = (set(), pd.DataFrame(columns=['Gene1', 'Gene2', 'y_true']))
    
    return results

def validate_network_files(datadir: Union[str, Path], datasets: List[Dict]) -> bool:
    """
    Validate that all specified network files exist and have correct format.
    
    :param datadir: Base data directory
    :param datasets: List of dataset configurations
    :returns: True if all files are valid
    """
    all_valid = True
    
    for dataset in datasets:
        dataset_name = dataset['name']
        true_edges = dataset.get('trueEdges')
        
        if isinstance(true_edges, str):
            network_files = [true_edges]
        elif isinstance(true_edges, list):
            network_files = true_edges
        else:
            print(f"❌ Invalid trueEdges format for dataset {dataset_name}: must be string or list")
            all_valid = False
            continue
        
        for network_file in network_files:
            network_path = Path(datadir) / dataset_name / network_file
            
            if not network_path.exists():
                print(f"❌ Network file not found: {network_path}")
                all_valid = False
                continue
            
            try:
                # Basic format validation
                df = pd.read_csv(network_path, sep=',', header=0, nrows=5)
                required_columns = ['Gene1', 'Gene2']
                
                if not all(col in df.columns for col in required_columns):
                    print(f"❌ Network file {network_path} missing required columns: {required_columns}")
                    all_valid = False
                    
            except Exception as e:
                print(f"❌ Failed to validate network file {network_path}: {e}")
                all_valid = False
    
    return all_valid 