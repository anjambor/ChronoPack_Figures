import os
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
from .singularityUtils import docker_to_singularity_cmd, get_singularity_image_path

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for GRNBoost2.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    # Define the base data directory relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    base_data_dir = project_root / 'BEELINE-data'

    # Construct the true input path by replacing the CWD with the BEELINE-data path
    # This is necessary because RunnerObj.inputDir incorrectly points to the CWD
    relative_input_path = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    true_input_dir = base_data_dir / relative_input_path

    # Ensure the output directory for GRNBOOST2 exists
    grnboost2_output_dir = RunnerObj.inputDir / "GRNBOOST2"
    grnboost2_output_dir.mkdir(parents=True, exist_ok=True)

    # Check if the final output file already exists
    if not (grnboost2_output_dir / "ExpressionData.csv").exists():
        # Construct the full path to the source expression data file
        expression_data_path = true_input_dir / RunnerObj.exprData
        
        # Read the expression data
        ExpressionData = pd.read_csv(expression_data_path, header=0, index_col=0)

        # Write the transposed data to the GRNBOOST2 output directory
        ExpressionData.T.to_csv(grnboost2_output_dir / "ExpressionData.csv", sep='\t', header=True, index=True)

    
def run(RunnerObj):
    '''
    Function to run GRNBOOST2 algorithm

    :param RunnerObj: An instance of the :class:`BLRun`
    '''
    # Get the BEELINE project root directory
    project_root = Path(__file__).resolve().parent.parent.parent

    # Construct paths relative to the project root
    # RunnerObj.inputDir is an absolute path, so make it relative to the project root
    relative_input_dir = RunnerObj.inputDir.relative_to(project_root)
    
    # The output directory is derived from the input directory's path relative to 'inputs'
    relative_dataset_path = relative_input_dir.relative_to(Path('inputs'))
    relative_output_dir = Path('outputs') / relative_dataset_path

    out_dir_host = project_root / relative_output_dir / "GRNBOOST2"
    os.makedirs(out_dir_host, exist_ok=True)

    # These paths are for inside the container and must be absolute
    in_file_container = f"/data/{relative_input_dir}/GRNBOOST2/ExpressionData.csv"
    out_file_container = f"/data/{relative_output_dir}/GRNBOOST2/outFile.txt"
    time_file_container = f"/data/{relative_output_dir}/GRNBOOST2/time.txt"

    # Get Singularity image path
    sif_path = get_singularity_image_path("GRNBOOST2")
    
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
        f'"DASK_TEMPORARY_DIRECTORY=/tmp '
        f'time -v -o {time_file_container} '
        f'python /runArboreto.py --algo=GRNBoost2 '
        f'--inFile={in_file_container} '
        f'--outFile={out_file_container}"'
    ]
    cmdToRun = ' '.join(cmd_parts)
    
    print(f"Running Singularity command: {cmdToRun}")
    return subprocess.Popen(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from GRNBOOST2.
    '''
    # Correctly determine the output directory relative to the project root
    project_root = Path(__file__).resolve().parent.parent.parent
    relative_input_dir = Path(str(RunnerObj.inputDir).replace(str(project_root), '').lstrip('/'))
    relative_dataset_path = relative_input_dir.relative_to('inputs')
    outDir = project_root / 'outputs' / relative_dataset_path / "GRNBOOST2"

    output_file = outDir / 'outFile.txt'

    if not output_file.exists():
        print(f"{output_file} does not exist, skipping...")
        return

    # Read the output data
    OutDF = pd.read_csv(output_file, sep='\t', header=0)
    
    # Write the ranked edges to a new CSV file
    ranked_edges_file = outDir / 'rankedEdges.csv'
    with open(ranked_edges_file, 'w') as outFile:
        outFile.write('Gene1\tGene2\tEdgeWeight\n')
        for idx, row in OutDF.iterrows():
            outFile.write(f"{row['TF']}\t{row['target']}\t{row['importance']}\n")
    
