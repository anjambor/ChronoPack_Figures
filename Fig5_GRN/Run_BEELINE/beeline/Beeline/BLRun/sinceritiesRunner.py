import os
import pandas as pd
from pathlib import Path
from .singularityUtils import get_singularity_image_path
import numpy as np
import subprocess


def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SINCERITIES.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Define the base data directory relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    base_data_dir = project_root / 'BEELINE-data'

    # Construct the true input path by replacing the CWD with the BEELINE-data path
    relative_input_path = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    true_input_dir = base_data_dir / relative_input_path

    # Ensure the output directory for SINCERITIES exists
    sincerities_output_dir = RunnerObj.inputDir / "SINCERITIES"
    sincerities_output_dir.mkdir(parents=True, exist_ok=True)

    # Construct full paths to the source data files
    expression_data_path = true_input_dir / RunnerObj.exprData
    cell_data_path = true_input_dir / RunnerObj.cellData

    # Read the expression and cell data, ensuring index_col=0 is used for consistency
    ExpressionData = pd.read_csv(expression_data_path, header=0, index_col=0)
    PTData = pd.read_csv(cell_data_path, header=0, index_col=0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = f"SINCERITIES/ExpressionData{idx}.csv"
        
        # Check if the output file already exists before processing
        if not (RunnerObj.inputDir / exprName).exists():
            newExpressionData = ExpressionData.loc[:,index].T
            
            # Perform quantile binning as recommended in the paper
            nBins = int(RunnerObj.params['nBins'])
            tQuantiles = pd.qcut(PTData.loc[index, colName], q=nBins, duplicates='drop')
            mid = [(a.left + a.right) / 2 for a in tQuantiles]

            newExpressionData['Time'] = mid
            newExpressionData.to_csv(RunnerObj.inputDir / exprName, sep=',', header=True, index=False)
    
def run(RunnerObj):
    '''
    Function to run SINCERITIES algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Get the BEELINE project root directory
    project_root = Path(__file__).resolve().parent.parent.parent

    # Construct paths relative to the project root
    relative_input_dir = RunnerObj.inputDir.relative_to(project_root)
    relative_dataset_path = relative_input_dir.relative_to(Path('inputs'))
    relative_output_dir = Path('outputs') / relative_dataset_path

    # Get Singularity image path
    sif_path = get_singularity_image_path("SINCERITIES")

    # Determine the number of time points by consistently reading cellData with index_col=0
    base_data_dir = project_root / 'BEELINE-data'
    true_input_dir = base_data_dir / relative_input_dir
    cell_data_path = true_input_dir / RunnerObj.cellData
    
    try:
        pt_data = pd.read_csv(cell_data_path, header=0, index_col=0)
        num_time_points = len(pt_data.columns)
    except FileNotFoundError:
        print(f"Warning: cellData file not found at {cell_data_path}. Assuming 1 time point.")
        num_time_points = 1

    # SINCERITIES handles multiple time points, so we loop through them
    for i in range(num_time_points):
        out_dir_host = project_root / relative_output_dir / "SINCERITIES"
        os.makedirs(out_dir_host, exist_ok=True)

        # These paths are for inside the container and must be absolute
        in_file_container = f"/data/{relative_input_dir}/SINCERITIES/ExpressionData{i}.csv"
        out_file_container = f"/data/{relative_output_dir}/SINCERITIES/outFile{i}.txt"
        time_file_container = f"/data/{relative_output_dir}/SINCERITIES/time{i}.txt"

        # Build Singularity command
        import tempfile
        temp_dir = tempfile.mkdtemp()
        
        cmd_parts = [
            'singularity', 'exec',
            '--cleanenv',
            '--bind', f'{project_root}:/data',
            '--bind', f'{temp_dir}:/tmp',
            str(sif_path),
            '/bin/sh', '-c',
            f'"cd /SINCERITIES && '
            f'time -v -o {time_file_container} '
            f'Rscript MAIN.R '
            f'{in_file_container} '
            f'{out_file_container}"'
        ]
        cmdToRun = ' '.join(cmd_parts)
        
        print(f"Running Singularity command: {cmdToRun}")
        p = subprocess.Popen(cmdToRun, shell=True)
        p.wait()

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SINCERITIES.

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Correctly determine the output and input directories relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    base_data_dir = project_root / 'BEELINE-data'
    
    relative_input_dir = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    true_input_dir = base_data_dir / relative_input_dir
    
    relative_dataset_path = relative_input_dir.relative_to('inputs')
    outDir = project_root / 'outputs' / relative_dataset_path / "SINCERITIES"

    # Consistently read cellData with index_col=0
    cell_data_path = true_input_dir / RunnerObj.cellData
    PTData = pd.read_csv(cell_data_path, header=0, index_col=0)
    colNames = PTData.columns
    OutSubDF = [0] * len(colNames)

    all_files_exist = True
    for idx in range(len(colNames)):
        outFile = f'outFile{idx}.txt'
        output_file = outDir / outFile
        if not output_file.exists():
            print(f"{output_file} does not exist, skipping...")
            all_files_exist = False
            break
        OutSubDF[idx] = pd.read_csv(output_file, sep=',', header=0)

    if not all_files_exist:
        return

    # Merge the dataframes by taking the maximum value from each DF
    outDF = pd.concat(OutSubDF)
    res = outDF[outDF['Interaction'] == outDF.groupby(['SourceGENES', 'TargetGENES'])['Interaction'].transform('max')]
    
    # Sort values and format the dataframe
    finalDF = res.sort_values('Interaction', ascending=False)
    finalDF.drop(labels='Edges', axis='columns', inplace=True)
    finalDF.columns = ['Gene2', 'Gene1', 'EdgeWeight']  # SINCERITIES output is incorrectly ordered
    
    # Save the ranked edges
    ranked_edges_file = outDir / 'rankedEdges.csv'
    finalDF.to_csv(ranked_edges_file, sep='\t', columns=['Gene1', 'Gene2', 'EdgeWeight'], index=False)
    
