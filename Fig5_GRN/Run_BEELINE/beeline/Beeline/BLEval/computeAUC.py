"""
Optimized AUPRC computation using Davis-Goadrich method.

This module now uses R's PRROC package for AUPRC calculation to ensure
consistency with computeDGAUC.py. The Davis-Goadrich method handles
interpolation and ties in precision-recall curves differently than
sklearn's standard AUC calculation.
"""

import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={"lines.linewidth": 2}, palette  = "deep", style = "ticks")
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from tqdm import tqdm
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector

# Import optimized utility functions
from .utils import get_truth_edges, get_pred_edges, make_truth_pred_df, profile_function, resolve_output_path

@profile_function
def PRROC(dataDict, inputSettings, directed = True, selfEdges = False, plotFlag = False, full_edge_universe = True, TFEdges = True):
    '''
    Computes areas under the precision-recall and ROC curves
    for a given dataset for each algorithm.
    
    OPTIMIZED VERSION: Uses efficient set operations and vectorized pandas
    operations instead of O(n⁴) dictionary-based approach.
    
    :param directed:   A flag to indicate whether to treat predictions as directed edges (directed = True) or undirected edges (directed = False).
    :type directed: bool
        
    :param selfEdges:   A flag to indicate whether to include self-edges (selfEdges = True) or exclude self-edges (selfEdges = False) from evaluation.
    :type selfEdges: bool
    
    :param plotFlag:   A flag to indicate whether or not to save PR and ROC plots.
    :type plotFlag: bool
    
    :param full_edge_universe:   A flag to control evaluation approach. If True, evaluates ALL possible edges like computeDGAUC.py (slower but matches exactly). If False, uses optimized approach with top-k edges only (faster).
    :type full_edge_universe: bool
    
    :param TFEdges:   A flag to indicate whether to consider only TF->gene edges (TFEdges = True) or all edges (TFEdges = False) for evaluation.
    :type TFEdges: bool
        
    :returns:
            - AUPRC: A dictionary containing AUPRC values for each algorithm
            - AUROC: A dictionary containing AUROC values for each algorithm
    '''
    
    # Load truth edges once using optimized function
    truth_path = str(inputSettings.datadir) + '/' + dataDict['name'] + '/' + dataDict['trueEdges']
    
    if full_edge_universe:
        # Load truth edges for full universe approach (similar to computeDGAUC.py)
        trueEdgesDF = pd.read_csv(truth_path, sep=',', header=0, index_col=None)
        
        # Apply TF filtering if requested
        if TFEdges:
            # Get TF genes (genes that appear as Gene1 in truth network)
            tf_genes = set(trueEdgesDF['Gene1'].unique())
            print(f"🧬 TF genes identified: {len(tf_genes)}")
            
            # Filter truth edges to only TF->gene edges
            trueEdgesDF = trueEdgesDF[trueEdgesDF['Gene1'].isin(tf_genes)]
            print(f"📊 Truth network (TF->gene only): {len(trueEdgesDF)} edges")
        
        unique_genes = np.unique(trueEdgesDF[['Gene1', 'Gene2']].values)
        print(f"📊 Truth network: {len(trueEdgesDF)} edges, {len(unique_genes)} unique genes")
        print(f"🌐 Full edge universe mode: Will evaluate ALL {len(unique_genes)}² possible edges")
    else:
        # Use optimized approach
        truth_set, truth_df = get_truth_edges(truth_path, directed=directed)
        
        # Apply TF filtering if requested
        if TFEdges:
            # Get TF genes (genes that appear as Gene1 in truth network)
            tf_genes = set(truth_df['Gene1'].unique())
            print(f"🧬 TF genes identified: {len(tf_genes)}")
            
            # Filter truth set to only TF->gene edges
            truth_set_tf = {(g1, g2) for g1, g2 in truth_set if g1 in tf_genes}
            truth_set = truth_set_tf
            
            # Filter truth DataFrame to only TF->gene edges
            truth_df = truth_df[truth_df['Gene1'].isin(tf_genes)]
            print(f"📊 Truth network (TF->gene only): {len(truth_set)} edges")
        
        print(f"📊 Truth network: {len(truth_set)} edges, {len(np.unique(truth_df[['Gene1', 'Gene2']].values))} unique genes")
        print(f"⚡ Optimized mode: Will evaluate only top-k predicted edges + truth edges")
            
    # Initialize result dictionaries
    precisionDict = {}
    recallDict = {}
    FPRDict = {}
    TPRDict = {}
    AUPRC = {}
    AUROC = {}
    
    # Process each algorithm
    for algo in tqdm(inputSettings.algorithms, 
                     total = len(inputSettings.algorithms), unit = " Algorithms"):

        algo_name = algo[0]
        
        # Skip algorithms that shouldn't run
        if not algo[1].get('should_run', True):
            print(f"⏸️  Skipping {algo_name} (should_run=False)")
            continue
        
        # Use smart path resolution (same as Early Precision)
        outDir = resolve_output_path(
            "outputs",  # Use hardcoded since we know this works
            str(inputSettings.datadir),
            dataDict['name'], 
            algo_name
        )
        pred_path = Path(outDir + '/rankedEdges.csv')
        
        # Check if the output rankedEdges file exists
        if not pred_path.exists():
            print(f"⚠️  {pred_path} does not exist. Skipping {algo_name}...")
            continue
        
        print(f"🔍 Processing {algo_name} from {pred_path}...")
        
        if full_edge_universe:
            # Full edge universe approach (matches computeDGAUC.py exactly)
            predEdgeDF = pd.read_csv(pred_path, sep='\t', header=0, index_col=None)
            
            if len(predEdgeDF) == 0:
                print(f"⚠️  No predictions found for {algo_name}")
                continue
            
            # Apply TF filtering to predictions if requested
            if TFEdges:
                # Filter predictions to only TF->gene edges
                predEdgeDF = predEdgeDF[predEdgeDF['Gene1'].isin(tf_genes)]
                print(f"📈 {algo_name}: {len(predEdgeDF)} TF->gene predicted edges loaded")
            else:
                print(f"📈 {algo_name}: {len(predEdgeDF)} predicted edges loaded")
            
            # Use the same logic as computeDGAUC.py
            eval_df = _build_full_edge_universe(trueEdgesDF, predEdgeDF, unique_genes, 
                                               directed=directed, selfEdges=selfEdges)
            
        else:
            # Optimized approach (current method)
            # Load predicted edges using optimized function
            # Limit to reasonable number of edges for evaluation efficiency
            
            # max_edges = min(len(truth_set) * 10, 100000)  # At most 10x truth edges or 100k
            max_edges = len(truth_set)
            pred_df = get_pred_edges(pred_path, top_k=max_edges, directed=directed)
            
            if len(pred_df) == 0:
                print(f"⚠️  No predictions found for {algo_name}")
                continue
            
            # Apply TF filtering to predictions if requested
            if TFEdges:
                # Filter predictions to only TF->gene edges
                pred_df = pred_df[pred_df['Gene1'].isin(tf_genes)]
                print(f"📈 {algo_name}: {len(pred_df)} TF->gene predicted edges loaded")
            else:
                print(f"📈 {algo_name}: {len(pred_df)} predicted edges loaded")
            
            # Create combined truth-prediction DataFrame for sklearn
            eval_df = make_truth_pred_df(truth_df, pred_df)
        
        print(f"🎯 {algo_name}: Evaluating {len(eval_df)} total edges ({eval_df['y_true'].sum()} true positives)")
        
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_true=eval_df['y_true'], 
                                y_score=eval_df['y_score'], 
                                pos_label=1)
        
        # Compute Precision-Recall curve  
        prec, recall, _ = precision_recall_curve(y_true=eval_df['y_true'],
                                                probas_pred=eval_df['y_score'], 
                                                pos_label=1)
        
        # Use R's PRROC package for Davis-Goadrich AUPRC calculation (same as computeDGAUC.py)
        try:
            prroc = importr('PRROC')
            prCurve = prroc.pr_curve(scores_class0 = FloatVector(list(eval_df['y_score'].values)), 
                      weights_class0 = FloatVector(list(eval_df['y_true'].values)))
            auprc_dg = prCurve[2][0]  # Davis-Goadrich AUPRC
        except Exception as e:
            print(f"⚠️  Warning: R PRROC not available for {algo_name}, falling back to sklearn: {e}")
            auprc_dg = auc(recall, prec)
        
        # Store results
        precisionDict[algo_name] = prec
        recallDict[algo_name] = recall
        FPRDict[algo_name] = fpr
        TPRDict[algo_name] = tpr
        AUPRC[algo_name] = auprc_dg  # Use Davis-Goadrich method
        AUROC[algo_name] = auc(fpr, tpr)
        
        print(f"✅ {algo_name}: AUPRC={AUPRC[algo_name]:.3f}, AUROC={AUROC[algo_name]:.3f}")
    
    # Generate plots if requested
    if plotFlag:
        _generate_plots(precisionDict, recallDict, FPRDict, TPRDict, 
                       AUPRC, AUROC, outDir, directed)
    
    print(f"🎯 PRROC completed: AUPRC has {len(AUPRC)} algorithms, AUROC has {len(AUROC)} algorithms")
    
    return AUPRC, AUROC

def _generate_plots(precisionDict, recallDict, FPRDict, TPRDict, 
                   AUPRC, AUROC, outDir, directed):
    """Generate PR and ROC plots for visualization."""
    
    PRName = '/PRplot' if directed else '/uPRplot'
    ROCName = '/ROCplot' if directed else '/uROCplot'
    
    # Make PR curves
    plt.figure(figsize=(10, 8))
    legendList = []
    for key in recallDict.keys():
        sns.lineplot(x=recallDict[key], y=precisionDict[key])
        legendList.append(key + ' (AUPRC = ' + str("%.3f" % (AUPRC[key])) + ')')
    
    plt.xlim(0, 1)    
    plt.ylim(0, 1)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curves')
    plt.legend(legendList, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(outDir + PRName + '.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(outDir + PRName + '.png', dpi=300, bbox_inches='tight')
    plt.clf()

    # Make ROC curves
    plt.figure(figsize=(10, 8))
    legendList = []
    for key in FPRDict.keys():
        sns.lineplot(x=FPRDict[key], y=TPRDict[key])
        legendList.append(key + ' (AUROC = ' + str("%.3f" % (AUROC[key])) + ')')

    plt.plot([0, 1], [0, 1], linewidth=1.5, color='k', linestyle='--', alpha=0.7)

    plt.xlim(0, 1)    
    plt.ylim(0, 1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves')
    plt.legend(legendList, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(outDir + ROCName + '.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(outDir + ROCName + '.png', dpi=300, bbox_inches='tight')
    plt.clf()

def _build_full_edge_universe(trueEdgesDF, predEdgeDF, unique_genes, directed=True, selfEdges=False):
    """
    Build the full edge universe approach identical to computeDGAUC.py.
    This creates ALL possible edges between genes and assigns truth/prediction values.
    
    :param trueEdgesDF: DataFrame with true edges
    :param predEdgeDF: DataFrame with predicted edges  
    :param unique_genes: Array of unique gene names
    :param directed: Whether to treat as directed graph
    :param selfEdges: Whether to include self-edges
    :returns: DataFrame with y_true and y_score columns for sklearn
    """
    from itertools import product, permutations, combinations, combinations_with_replacement
    
    if directed:
        # Initialize dictionaries with all possible edges
        if selfEdges:
            possibleEdges = list(product(unique_genes, repeat=2))
        else:
            possibleEdges = list(permutations(unique_genes, r=2))
        
        TrueEdgeDict = {'|'.join(p): 0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p): 0 for p in possibleEdges}
        
        # Compute TrueEdgeDict Dictionary
        # 1 if edge is present in the ground-truth, 0 otherwise
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[(trueEdgesDF['Gene1'] == key.split('|')[0]) &
                   (trueEdgesDF['Gene2'] == key.split('|')[1])]) > 0:
                TrueEdgeDict[key] = 1
                
        # Compute PredEdgeDict Dictionary  
        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[(predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])]
            if len(subDF) > 0:
                PredEdgeDict[key] = np.abs(subDF.EdgeWeight.values[0])

    else:
        # For undirected graphs
        if selfEdges:
            possibleEdges = list(combinations_with_replacement(unique_genes, r=2))
        else:
            possibleEdges = list(combinations(unique_genes, r=2))
            
        TrueEdgeDict = {'|'.join(p): 0 for p in possibleEdges}
        PredEdgeDict = {'|'.join(p): 0 for p in possibleEdges}

        # Compute TrueEdgeDict Dictionary
        for key in TrueEdgeDict.keys():
            if len(trueEdgesDF.loc[((trueEdgesDF['Gene1'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene2'] == key.split('|')[1])) |
                              ((trueEdgesDF['Gene2'] == key.split('|')[0]) &
                           (trueEdgesDF['Gene1'] == key.split('|')[1]))]) > 0:
                TrueEdgeDict[key] = 1  

        # Compute PredEdgeDict Dictionary
        for key in PredEdgeDict.keys():
            subDF = predEdgeDF.loc[((predEdgeDF['Gene1'] == key.split('|')[0]) &
                               (predEdgeDF['Gene2'] == key.split('|')[1])) |
                              ((predEdgeDF['Gene2'] == key.split('|')[0]) &
                               (predEdgeDF['Gene1'] == key.split('|')[1]))]
            if len(subDF) > 0:
                PredEdgeDict[key] = max(np.abs(subDF.EdgeWeight.values))
                
    # Combine into one dataframe to pass to sklearn/PRROC
    outDF = pd.DataFrame([TrueEdgeDict, PredEdgeDict]).T
    outDF.columns = ['y_true', 'y_score']
    
    return outDF
