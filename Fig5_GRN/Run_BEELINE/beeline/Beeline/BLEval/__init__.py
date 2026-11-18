"""
BEELINE Evaluation (:mod:`BLEval`) module contains the following
:class:`BLEval.BLEval` and three additional classes used in the
definition of BLEval class 

- :class:`BLEval.ConfigParser` 
- :class:`BLEval.InputSettings` 
- :class:`BLEval.OutputSettings`


"""
import os
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from glob import glob
import pathlib
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency


# local imports
from BLEval.parseTime import getTime
# from BLEval.computeDGAUC import PRROC
from BLEval.computeAUC import PRROC
from BLEval.computeBorda import Borda
from BLEval.computeJaccard import Jaccard
from BLEval.computeSpearman import Spearman
from BLEval.computeNetMotifs import Motifs
from BLEval.computeEarlyPrec import EarlyPrec
from BLEval.computePathStats import pathAnalysis
from BLEval.computeSignedEPrec import signedEPrec


class InputSettings(object):
    '''
    The class for storing the names of input files.
    This initilizes an InputSettings object based on the
    following parameters.
    

    :param datadir:   input dataset root directory, typically 'inputs/'
    :type datadir: str

    :param datasets:   List of dataset configurations (original format)
    :type datasets: list

    :param network_configs:   List of expanded network configurations for evaluation
    :type network_configs: list
        
    :param algorithms:   List of algorithm names
    :type algorithms: list
    '''

    def __init__(self,
            datadir, datasets, network_configs, algorithms) -> None:

        self.datadir = datadir
        self.datasets = datasets  # Original dataset configs for backward compatibility
        self.network_configs = network_configs  # Expanded network configs for evaluation
        self.algorithms = algorithms


class OutputSettings(object):
    '''    
    The class for storing the names of directories that output should
    be written to. This initilizes an OutputSettings object based on the
    following two parameters.
    

    :param base_dir: output root directory, typically 'outputs/'
    :type base_dir: str
    :param output_prefix: A prefix added to the final output files.
    :type str:
    '''

    def __init__(self, base_dir, output_prefix: Path) -> None:
        # Smart path resolution: check if the directory exists before prepending "../"
        base_dir_path = Path(base_dir)
        
        if not base_dir_path.is_absolute():
            # First try the path as-is (for when running from project root)
            if base_dir_path.exists():
                self.base_dir = str(base_dir_path)
            # If not found and doesn't start with "../", try with parent reference
            elif not str(base_dir_path).startswith("../"):
                parent_path = Path("..").joinpath(base_dir_path)
                if parent_path.exists():
                    self.base_dir = str(parent_path)
                else:
                    # Fall back to original path if neither exists
                    self.base_dir = str(base_dir_path)
            else:
                self.base_dir = str(base_dir_path)
        else:
            self.base_dir = str(base_dir_path)
            
        self.output_prefix = output_prefix




class BLEval(object):
    '''
    The BEELINE Evaluation object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs.
    '''

    def __init__(self,
            input_settings: InputSettings,
            output_settings: OutputSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings


    def computeAUC(self, directed = True, selfEdges = False, full_edge_universe = False, TFEdges = False):

        '''
        Computes areas under the precision-recall (PR) and
        and ROC plots for each algorithm-network combination.
        
        Includes a "Random" baseline row representing the expected performance
        of a random predictor for comparison.
      
        Parameters
        ----------
        directedFlag: bool
            A flag to specifiy whether to treat predictions
            as directed edges (directed = True) or 
            undirected edges (directed = False).
        selfEdges: bool
            A flag to indicate whether to include self-edges 
            (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
        full_edge_universe: bool
            A flag to control evaluation approach. If True, evaluates ALL possible edges 
            like computeDGAUC.py (slower but matches exactly). If False, uses optimized 
            approach with top-k edges only (faster).
        TFEdges: bool
            A flag to indicate whether to consider only TF->gene edges (TFEdges = True) 
            or all edges (TFEdges = False) for evaluation.
        
        :returns:
            - AUPRC: A dataframe containing AUPRC values for each algorithm-dataset combination
            - AUROC: A dataframe containing AUROC values for each algorithm-dataset combination
        '''
        from .utils import get_truth_edges
        
        AUPRCDict = {}
        AUROCDict = {}

        # Use expanded network configurations for evaluation
        for network_config in tqdm(self.input_settings.network_configs, 
                                  total = len(self.input_settings.network_configs), 
                                  unit = " Networks"):
            
            # Create a dataset-like object for backward compatibility with PRROC function
            dataset_obj = {
                'name': network_config['dataset_name'],
                'trueEdges': network_config['network_file']
            }
            
            AUPRC, AUROC = PRROC(dataset_obj, self.input_settings, 
                                    directed = directed, selfEdges = selfEdges, plotFlag = False, 
                                    full_edge_universe = full_edge_universe, TFEdges = TFEdges)
            
            # Use column name for results (e.g., "K562_ChIP-seq")
            column_name = network_config['column_name']
            AUPRCDict[column_name] = AUPRC
            AUROCDict[column_name] = AUROC
            
            # Add random baseline for this network
            truth_path = str(self.input_settings.datadir) + '/' + network_config['dataset_name'] + '/' + network_config['network_file']
            truth_set, truth_df = get_truth_edges(truth_path, directed=directed)
            
            # Apply TF filtering if requested (same as in PRROC function)
            if TFEdges:
                tf_genes = set(truth_df['Gene1'].unique())
                truth_set_tf = {(g1, g2) for g1, g2 in truth_set if g1 in tf_genes}
                truth_set = truth_set_tf
                truth_df = truth_df[truth_df['Gene1'].isin(tf_genes)]
            
            unique_genes = np.unique(truth_df[['Gene1', 'Gene2']].values)
            
            # Calculate network density (prevalence of true edges)
            if TFEdges:
                # For TF->gene edges, only TF genes can be sources
                tf_genes = set(truth_df['Gene1'].unique()) if len(truth_df) > 0 else set()
                target_genes = set(truth_df['Gene2'].unique()) if len(truth_df) > 0 else set()
                
                if directed:
                    if selfEdges:
                        # TF can regulate itself: |TF| * |all_genes|
                        n_possible = len(tf_genes) * len(unique_genes)
                    else:
                        # TF cannot regulate itself: |TF| * (|all_genes| - 1) for self-TFs, |TF| * |all_genes| for non-self
                        n_possible = len(tf_genes) * len(unique_genes) - len(tf_genes.intersection(unique_genes))
                else:
                    # Undirected TF edges don't make biological sense, but handle it
                    if selfEdges:
                        n_possible = len(tf_genes) * (len(tf_genes) + 1) // 2
                    else:
                        n_possible = len(tf_genes) * (len(tf_genes) - 1) // 2
            else:
                # Standard calculation for all edges
                if directed:
                    if selfEdges:
                        n_possible = len(unique_genes) ** 2
                    else:
                        n_possible = len(unique_genes) * (len(unique_genes) - 1)
                else:
                    if selfEdges:
                        n_possible = len(unique_genes) * (len(unique_genes) + 1) // 2
                    else:
                        n_possible = len(unique_genes) * (len(unique_genes) - 1) // 2
            
            network_density = len(truth_set) / n_possible if n_possible > 0 else 0.0
            
            # Random predictor baseline
            # AUROC: always 0.5 for random predictor (independent of class imbalance)
            # AUPRC: depends on evaluation approach
            if full_edge_universe:
                # Full universe: prevalence = network density
                baseline_auprc = network_density
            else:
                # Optimized mode: evaluation pool contains k true edges + ~k predicted edges
                # Expected prevalence ≈ k/(k+k) = 0.5 when algorithm gets no true edges right
                # This is the correct baseline for the actual evaluation pool used
                baseline_auprc = 0.5
            
            # Add to dictionaries - use algorithm name that won't conflict
            baseline_name = "Random"
            if baseline_name in AUPRCDict[column_name]:
                baseline_name = "RandomPredictor"
            
            AUPRCDict[column_name][baseline_name] = baseline_auprc
            AUROCDict[column_name][baseline_name] = 0.5
        
        # Convert to DataFrame - handle both single and multiple datasets
        print(f"🔄 Converting results to DataFrame: {len(AUPRCDict)} datasets")
        
        if AUPRCDict and any(AUPRCDict.values()):  # Check if we have any results
            AUPRC_df = pd.DataFrame(AUPRCDict)
            AUROC_df = pd.DataFrame(AUROCDict)
        else:
            # Create empty DataFrames with proper structure if no results
            AUPRC_df = pd.DataFrame()
            AUROC_df = pd.DataFrame()
            
        print(f"📊 DataFrame shapes - AUPRC: {AUPRC_df.shape}, AUROC: {AUROC_df.shape}")
        print(f"🎲 Random baseline added for all datasets")
        return AUPRC_df, AUROC_df
    

    def computeEarlyPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values of the 
        network formed using the predicted top-k edges.
        
        Includes a "Random" baseline row representing the expected performance
        of a random predictor for comparison.
        
        :returns: 
            A dataframe containing the early precision values 
            for each algorithm-dataset combination.

        '''
        from .utils import get_truth_edges
        
        Eprec = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                Eprec[algo[0]] = EarlyPrec(self, algo[0])
        
        # Add random baseline row
        # For early precision with random predictor: expected precision = network density
        baseline_row = {}
        for network_config in self.input_settings.network_configs:
            truth_path = str(self.input_settings.datadir) + '/' + network_config['dataset_name'] + '/' + network_config['network_file']
            truth_set, truth_df = get_truth_edges(truth_path, directed=True)
            unique_genes = np.unique(truth_df[['Gene1', 'Gene2']].values)
            
            # Calculate network density (directed, no self-edges)
            # Note: Early precision always uses directed=True, so we match that here
            n_possible = len(unique_genes) * (len(unique_genes) - 1)
            network_density = len(truth_set) / n_possible if n_possible > 0 else 0.0
            baseline_row[network_config['column_name']] = network_density
        
        # Use algorithm name that won't conflict
        baseline_name = "Random"
        if baseline_name in Eprec:
            baseline_name = "RandomPredictor"
        
        Eprec[baseline_name] = baseline_row
        
        print(f"🎲 Random baseline added to Early Precision results")
        return pd.DataFrame(Eprec).T


    def parseTime(self):
        """
        Parse time output for each
        algorithm-dataset combination.
        
        Note: Since algorithms run once per dataset (not per network),
        time results are still organized by dataset name only.
        
        :returns:
            A dictionary of times for all dataset-algorithm combinations
        """
        TimeDict = dict()

        # Use original datasets for timing since algorithms run per dataset
        for dataset in tqdm(self.input_settings.datasets, 
                            total = len(self.input_settings.datasets), unit = " Datasets"):
            timevals  = getTime(self, dataset)
            TimeDict[dataset["name"]] = timevals

        return TimeDict

    def computeJaccard(self):

        '''
        Computes Jaccard Index between top-k edge predictions
        of the same algorithm.
        
        :returns:
            A dataframe containing the median and median absolute
            deviation of the Jaccard Index values of each algorithm
            on the given set of datasets.
        '''
        JaccDF = {}
        JaccDF['Jaccard Median'] = {}
        JaccDF['Jaccard MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                JaccDF['Jaccard Median'][algo[0]], JaccDF['Jaccard MAD'][algo[0]] = Jaccard(self, algo[0])
            
        return pd.DataFrame(JaccDF)     

    
    def computeSpearman(self):

        '''
        Finds the Spearman's correlation coefficient between
        the ranked edges of the same algorithm on the given
        set of datasets.
        
        :returns:
            A dataframe containing the median and median absolute
            deviation of the Separman's correlation values of each algorithm.
        '''
        corrDF = {}
        corrDF['Spearman Median'] = {}
        corrDF['Spearman MAD'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                corrDF['Spearman Median'][algo[0]],corrDF['Spearman MAD'][algo[0]] = Spearman(self, algo[0])
            
        return pd.DataFrame(corrDF)



    def computeNetMotifs(self):

        '''
        For each algorithm-network combination, this function computes the network motifs such as 
        Feedforward loops, Feedback loops and 
        Mutual interactions in the predicted top-k network. It returns the ratio of network motif counts
        compared to their respective values in the reference network.
      
        :returns:
            - FBL: A dataframe containing ratios of number of Feedback loops
            
            - FFL: A dataframe containing ratios of number of Feedforward loops
            
            - MI: A dataframe containing ratios of number of Mutual Interactions
        '''
        FFLDict = {}
        FBLDict = {}
        MIDict = {}
        
        # Use expanded network configurations
        for network_config in tqdm(self.input_settings.network_configs, 
                                  total = len(self.input_settings.network_configs), unit = " Networks"):
            # Create a dataset-like object for backward compatibility
            dataset_obj = {
                'name': network_config['dataset_name'],
                'trueEdges': network_config['network_file']
            }
            column_name = network_config['column_name']
            
            FBLDict[column_name], FFLDict[column_name], MIDict[column_name] = Motifs(dataset_obj, self.input_settings)
            
        FBL = pd.DataFrame(FBLDict)
        FFL = pd.DataFrame(FFLDict)
        MI = pd.DataFrame(MIDict)
        
        return FBL, FFL, MI
    
    
    def computePaths(self):
        '''
        For each algorithm-network combination, this function computes path lengths
        through TP edges and FP edges, returns statistics on path lengths.

        :returns:
            - pathStats: A dataframe path lengths in predicted network
        '''
        # Use expanded network configurations
        for network_config in tqdm(self.input_settings.network_configs, 
                                  total = len(self.input_settings.network_configs), unit = " Networks"):
            # Create a dataset-like object for backward compatibility
            dataset_obj = {
                'name': network_config['dataset_name'],
                'trueEdges': network_config['network_file']
            }
            pathAnalysis(dataset_obj, self.input_settings)
    
    def computeSignedEPrec(self):

        '''
        For each algorithm-dataset combination,
        this function computes the Early Precision values separately 
        for the activation and inhibitory edges.
        
        :returns: 
            - A dataframe containing early precision for activation edges
            - A dataframe containing early precision for inhibitory edges
        '''
        sEPRDict = {}
        sEPRDict['EPrec Activation'] = {}
        sEPRDict['EPrec Inhibition'] = {}
        outDir = str(self.output_settings.base_dir) + \
                 str(self.input_settings.datadir).split("inputs")[1] + "/"
        for algo in tqdm(self.input_settings.algorithms, unit = " Algorithms"):
            if algo[1]['should_run'] == True:
                sEPrecDF = signedEPrec(self, algo[0])
                sEPRDict['EPrec Activation'][algo[0]] = sEPrecDF['+']
                sEPRDict['EPrec Inhibition'][algo[0]] = sEPrecDF['-']
        return(pd.DataFrame(sEPRDict['EPrec Activation']).T, pd.DataFrame(sEPRDict['EPrec Inhibition']).T)
    
    def computeBorda(self, selectedAlgorithms=None, aggregationMethod="average"):

        '''
        Computes edge ranked list using the Borda method for each dataset.
        Parameters
        ----------
        selectedAlgorithms: [str]
          List of algorithm names used to run borda method on selected
          algorithms. If nothing is provided, the function runs borda on
          all available algorithms.
        aggregationMethod: str
          Method used to aggregate rank in borda method. Available options are
          {‘average’, ‘min’, ‘max’, ‘first’}, default ‘average’
        :returns:
            None
        '''
        feasibleAlgorithmOptions = [algorithmName for algorithmName, _ in self.input_settings.algorithms]
        feasibleaggregationMethodOptions = ['average', 'min', 'max', 'first']

        selectedAlgorithms = feasibleAlgorithmOptions if selectedAlgorithms is None else selectedAlgorithms

        for a in selectedAlgorithms:
            if a not in feasibleAlgorithmOptions:
                print("\nERROR: No data available on algorithm %s. Please choose an algorithm from the following options: %s" % (a, feasibleAlgorithmOptions))
                return

        if aggregationMethod not in feasibleaggregationMethodOptions:
                print("\nERROR: Please choose an aggregation method algorithm from following options: " % feasibleaggregationMethodOptions)
                return

        Borda(self, selectedAlgorithms, aggregationMethod)

class ConfigParser(object):
    '''
    The class define static methods for parsing and storing the contents 
    of the config file that sets a that sets a large number of parameters 
    used in the BLEval.
    '''
    @staticmethod
    def parse(config_file_handle) -> BLEval:
        '''
        A method for parsing the input .yaml file.
        
        :param config_file_handle: Name of the .yaml file to be parsed
        :type config_file_handle: str
        
        :returns: 
            An object of class :class:`BLEval.BLEval`.

        '''
        config_map = yaml.load(config_file_handle)
        return BLEval(
            ConfigParser.__parse_input_settings(
                config_map['input_settings']),
            ConfigParser.__parse_output_settings(
                config_map['output_settings']))
    
    @staticmethod
    def __parse_input_settings(input_settings_map) -> InputSettings:
        '''
        A method for parsing and initializing 
        InputSettings object with support for multiple ground truth networks.
        '''
        from .utils import parse_multiple_networks, validate_network_files
        
        input_dir = input_settings_map['input_dir']
        dataset_dir = input_settings_map['dataset_dir']
        
        # Convert datadir to absolute path rooted at project directory
        datadir_path = Path(input_dir, dataset_dir)
        if not datadir_path.is_absolute():
            # Get BEELINE project root (parent of the Beeline subdirectory) - same as BLRunner
            project_root = Path(__file__).resolve().parent.parent.parent
            # The actual data is in BEELINE-data/inputs, not project root inputs
            datadir_path = (project_root / "BEELINE-data" / datadir_path).resolve()
        
        # Check if datasets have been specified or not
        if 'datasets' in input_settings_map: 
            datasets_specified = True
        else:
            datasets_specified = False

        # If no datasets specified, run all datasets in dataset_dir
        if datasets_specified is False: 
             subfolder_dir = glob(os.path.join(str(datadir_path), "*/"), recursive = True)
             datasets = []
             for x in subfolder_dir:
                datasets.append({"name": pathlib.Path(x).name, 
                                 "exprData": "ExpressionData.csv", 
                                 "cellData": "PseudoTime.csv", 
                                 "trueEdges": "refNetwork.csv"})
        # If datasets specified, run the corresponding datasets
        else:
            datasets = input_settings_map['datasets']
            if datasets is None:
                print("Please specify input datasets!")

        # Validate all network files exist
        print("🔍 Validating ground truth network files...")
        if not validate_network_files(datadir_path, datasets):
            print("❌ Some network files are missing or invalid. Please check the configuration.")
        
        # Expand datasets into individual network configurations
        network_configs = []
        for dataset in datasets:
            try:
                configs = parse_multiple_networks(dataset)
                network_configs.extend(configs)
            except Exception as e:
                print(f"❌ Error parsing networks for dataset {dataset.get('name', 'unknown')}: {e}")
                
        print(f"📊 Expanded {len(datasets)} datasets into {len(network_configs)} network evaluation configurations")
        for config in network_configs:
            print(f"   • {config['column_name']} ({config['network_file']})")

        return InputSettings(
                datadir_path,
                datasets,
                network_configs,
                ConfigParser.__parse_algorithms(
                input_settings_map['algorithms']))


    @staticmethod
    def __parse_algorithms(algorithms_list):
        '''
        A method for parsing the list of algorithms
        that are being evaluated, along with
        any parameters being passed.
        
        Note that these parameters may not be
        used in the current evaluation, but can 
        be used at a later point.
        '''
        
        # Initilalize the list of algorithms
        algorithms = []
        
        # Parse contents of algorithms_list
        for algorithm in algorithms_list:
                combos = [dict(zip(algorithm['params'], val))
                    for val in itertools.product(
                        *(algorithm['params'][param]
                            for param in algorithm['params']))]
                for combo in combos:
                    algorithms.append([algorithm['name'],combo])
            

        return algorithms

    @staticmethod
    def __parse_output_settings(output_settings_map) -> OutputSettings:
        '''
        A method for parsing and initializing 
        Output object.
        '''
        output_dir = Path(output_settings_map['output_dir'])
        output_prefix = Path(output_settings_map['output_prefix'])

        return OutputSettings(output_dir,
                             output_prefix)
