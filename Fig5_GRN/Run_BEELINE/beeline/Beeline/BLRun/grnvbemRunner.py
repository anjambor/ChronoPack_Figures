import os
import pandas as pd
from pathlib import Path
from .singularityUtils import get_singularity_image_path
import numpy as np
import tempfile
import subprocess

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRNVBEM.
    It will create the input folder at RunnerObj.dataset.name/GRNVBEM/ if it
    does not exist already. The input folder will contain an ExpressionData.csv with
    cells ordered according to the pseudotime along the columns, and genes along
    the rows. If the files already exist, this function will overwrite it.
    
    '''
    # Define the base data directory relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    base_data_dir = project_root / 'BEELINE-data'

    # Construct the true input path by replacing the CWD with the BEELINE-data path
    relative_input_path = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    true_input_dir = base_data_dir / relative_input_path

    # Ensure the output directory for GRNVBEM exists
    grnvbem_output_dir = RunnerObj.inputDir / "GRNVBEM"
    grnvbem_output_dir.mkdir(parents=True, exist_ok=True)

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
        exprName = f"GRNVBEM/ExpressionData{idx}.csv"
        
        # Check if the output file already exists before processing
        if not (RunnerObj.inputDir / exprName).exists():
            subPT = PTData.loc[index,:]
            subExpr = ExpressionData[index]
            # Order columns by PseudoTime
            newExpressionData = subExpr[subPT.sort_values([colName]).index.astype(str)]
            
            newExpressionData.insert(loc=0, column='GENES', value=newExpressionData.index)

            # Write .csv file
            newExpressionData.to_csv(RunnerObj.inputDir / exprName, sep=',', header=True, index=False)

    
def run(RunnerObj):
    '''
    Function to run GRNVBEM algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Get the BEELINE project root directory
    project_root = Path(__file__).resolve().parent.parent.parent

    # Construct paths relative to the project root
    relative_input_dir = RunnerObj.inputDir.relative_to(project_root)
    relative_dataset_path = relative_input_dir.relative_to(Path('inputs'))
    relative_output_dir = Path('outputs') / relative_dataset_path

    # Get Singularity image path
    sif_path = get_singularity_image_path("GRNVBEM")

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

    # GRNVBEM handles multiple time points, so we loop through them
    for i in range(num_time_points):
        out_dir_host = project_root / relative_output_dir / "GRNVBEM"
        os.makedirs(out_dir_host, exist_ok=True)

        # These paths are for inside the container and must be absolute
        in_file_container = f"/data/{relative_input_dir}/GRNVBEM/ExpressionData{i}.csv"
        out_file_container = f"/data/{relative_output_dir}/GRNVBEM/outFile{i}.txt"
        time_file_container = f"/data/{relative_output_dir}/GRNVBEM/time{i}.txt"

        # Build Singularity command
        temp_dir = tempfile.mkdtemp()
        
        cmd_parts = [
            'singularity', 'exec',
            '--cleanenv',
            '--bind', f'{project_root}:/data',
            '--bind', f'{temp_dir}:/tmp',
            str(sif_path),
            '/bin/sh', '-c',
            f'"time -v -o {time_file_container} '
            f'/VBEM/GRNVBEM '
            f'{in_file_container} '
            f'{out_file_container}"'
        ]
        cmdToRun = ' '.join(cmd_parts)
        
        print(f"Running Singularity command: {cmdToRun}")
        p = subprocess.Popen(cmdToRun, shell=True)
        p.wait()

def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRNVBEM.
    '''
    # Correctly determine the output and input directories relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    base_data_dir = project_root / 'BEELINE-data'
    
    relative_input_dir = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    true_input_dir = base_data_dir / relative_input_dir
    
    relative_dataset_path = relative_input_dir.relative_to('inputs')
    outDir = project_root / 'outputs' / relative_dataset_path / "GRNVBEM"

    # Consistently read cellData with index_col=0
    cell_data_path = true_input_dir / RunnerObj.cellData
    PTData = pd.read_csv(cell_data_path, header=0, index_col=0)

    colNames = PTData.columns
    OutSubDF = [0] * len(colNames)

    all_files_exist = True
    for indx in range(len(colNames)):
        outFileName = f'outFile{indx}.txt'
        output_file = outDir / outFileName
        if not output_file.exists():
            print(f"{output_file} does not exist, skipping...")
            all_files_exist = False
            break
        
        # Read output
        OutSubDF[indx] = pd.read_csv(output_file, sep='\t', header=0)
        
    if not all_files_exist:
        return

    outDF = pd.concat(OutSubDF)
    FinalDF = outDF[outDF['Probability'] == outDF.groupby(['Parent', 'Child'])['Probability'].transform('max')]

    ranked_edges_file = outDir / 'rankedEdges.csv'
    with open(ranked_edges_file, 'w') as outFile:
        outFile.write('Gene1\tGene2\tEdgeWeight\n')
        for idx, row in FinalDF.sort_values(['Probability'], ascending=False).iterrows():
            outFile.write(f"{row['Parent']}\t{row['Child']}\t{row['Probability']}\n")
    
