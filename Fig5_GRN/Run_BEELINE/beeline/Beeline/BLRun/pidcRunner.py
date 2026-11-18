import os
import pandas as pd
from pathlib import Path
import numpy as np
from .singularityUtils import get_singularity_image_path
import subprocess

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SCODE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("PIDC").exists():
        print("Input folder for PIDC does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("PIDC").mkdir(parents=True, exist_ok=True)
        
    if not RunnerObj.inputDir.joinpath("PIDC/ExpressionData.csv").exists():
        ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
        ExpressionData.to_csv(RunnerObj.inputDir.joinpath("PIDC/ExpressionData.csv"),
                             sep = '\t', header  = True, index = True)
    
def run(RunnerObj):
    '''
    Function to run PIDC algorithm using Singularity
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/PIDC/ExpressionData.csv"
    
    # make output dirs if they do not exist:
    outDir = "../outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PIDC/"
    os.makedirs(outDir, exist_ok = True)
    
    outPath = 'data/'+ str(outDir) + 'outFile.txt'
    
    # Get Singularity image path
    sif_path = get_singularity_image_path("PIDC")
    
    # Build Singularity command - ensure proper container execution with working directory
    # Add temp directory for Julia cache operations while keeping pre-installed packages accessible
    import tempfile
    temp_dir = tempfile.mkdtemp()
    cmdToRun = ' '.join([
        'singularity exec',
        '--bind', f'{str(Path.cwd())}:/data',
        '--bind', f'{temp_dir}:/tmp',
        '--pwd', '/',
        str(sif_path),
        '/bin/sh', '-c',
        f'"JULIA_DEPOT_PATH=/root/.julia:/tmp time -v -o data/{outDir}time.txt julia /runPIDC.jl {inputPath} {outPath}"'
    ])
    
    print(f"Running Singularity command: {cmdToRun}")
    return subprocess.Popen(cmdToRun, shell=True)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SCODE.
    '''
    # Quit if output directory does not exist
    outDir = "../outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/PIDC/"
    if not Path(outDir+'outFile.txt').exists():
        print(outDir+'outFile.txt'+'does not exist, skipping...')
        return
        
    # Read output
    OutDF = pd.read_csv(outDir+'outFile.txt', sep = '\t', header = None)
    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in OutDF.iterrows():
        outFile.write('\t'.join([row[0],row[1],str(row[2])])+'\n')
    outFile.close()
    
