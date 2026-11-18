import os
import subprocess
import pandas as pd
from pathlib import Path
from .singularityUtils import get_singularity_image_path
import numpy as np

def generateInputs(RunnerObj):
    '''
    Function to generate desired inputs for SINGE.
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("SINGE").exists():
        print("Input folder for SINGE does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("SINGE").mkdir(parents=True, exist_ok=True)
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "SINGE/ExpressionData"+str(idx)+".csv"
        newExpressionData = ExpressionData.loc[:,index].T
        newExpressionData['PseudoTime'] = PTData.loc[index,colName]
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)


def run(RunnerObj):
    '''
    Function to run SINGE algorithm
    '''
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1] + \
                    "/SINGE/"
    

    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINGE/"
    os.makedirs(outDir, exist_ok = True)

    # if the parameters aren't specified, then use default parameters
    # TODO allow passing in multiple sets of hyperparameters
    # these must be in the right order!
    params_order = [
        'lambda', 'dT', 'num_lags', 'kernel_width',
        'prob_zero_removal', 'prob_remove_samples',
        'family'
    ]
    default_params = {
        'lambda': '0.01',
        'dT': '10',
        'num_lags': '5',
        'kernel_width': '4',
        'prob_zero_removal': '0',
        'prob_remove_samples': '0.2',
        'family': 'gaussian',
        'num_replicates': '2',
    }
    params = RunnerObj.params
    for param, val in default_params.items():
        if param not in params:
            params[param] = val
    
    num_replicates = params['num_replicates']
    replicates = []
    for replicate in range(num_replicates):
       replicates.append(' '.join('--' + p.replace('_', '-') + ' ' + str(params[p]) for p in params_order) + ' '.join(['', '--replicate', str(replicate), '--ID', str(replicate)]))
    params_str = '\n'.join(replicates)    

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):    
        # Ensure SINGE input directory exists
        os.makedirs(str(RunnerObj.inputDir.joinpath("SINGE")), exist_ok=True)
        
        outPath = str(outDir) + str(idx) + "/" 
        os.makedirs(outPath, exist_ok = True)
        outFile = "data/" + outPath
        outFileSymlink = "out" + str(idx)
        inputFile = inputPath + "ExpressionData"+str(idx)+".csv"
        inputMat = inputPath + "ExpressionData"+str(idx)+".mat"
        geneListMat = inputPath + "GeneList"+str(idx)+".mat"
        paramsFile = inputPath + "hyperparameters.txt"       

        '''
        This is a workaround for https://github.com/gitter-lab/SINGE/blob/master/code/parseParams.m#L39
        not allowing '/' characters in the outDir parameter.
        '''
        symlink_out_file = ' '.join(['ln -s', '/usr/local/SINGE/' + outFile, outFileSymlink])

        '''
        See https://github.com/gitter-lab/SINGE/blob/master/README.md.  SINGE expects a data matfile with variables "X" and "ptime",
        and a gene_list matfile with the variable "gene_list".

        Saving fullKp is a very hacky workaround for https://github.com/gitter-lab/SINGE/blob/master/code/iLasso_for_SINGE.m#L56, 
        that assumes this input was saved in matfile v7.3 which octave does not support.
        '''
        # Create octave command with proper escaping
        octave_script = f"CSV = csvread('/usr/local/SINGE/{inputFile}'); " + \
                       "X = sparse(CSV(2:end,1:end-1).'); ptime = CSV(2:end,end).'; " + \
                       "Kp2.Kp = single(ptime); Kp2.sumKp = single(ptime*X.'); " + \
                       f"fullKp(1, {params['dT']*params['num_lags']}) = Kp2; " + \
                       f"save('-v7','/usr/local/SINGE/{inputMat}', 'X', 'ptime', 'fullKp'); " + \
                       f"f = fopen('/usr/local/SINGE/{inputFile}'); " + \
                       "gene_list = strsplit(fgetl(f), ','); gene_list = gene_list(1:end-1).'; fclose(f); " + \
                       f"save('-v7','/usr/local/SINGE/{geneListMat}', 'gene_list');"

        cmdToRun = ' '.join(['singularity exec --fakeroot --bind', 
                             str(Path.cwd())+':/usr/local/SINGE/data', 
                             str(get_singularity_image_path('grnbeeline/singe:0.4.1')),
                             '/bin/sh -c "echo \\"' + params_str + '\\" > /usr/local/SINGE/' + paramsFile,
                             '&&', symlink_out_file,
                             '&& octave -q --eval \\"' + octave_script + '\\"',
                             '&& mkdir -p /usr/local/SINGE/data/' + str(outDir),
                             '&& time -v -o /usr/local/SINGE/data/' + str(outDir) + 'time'+str(idx)+'.txt',
                             '/usr/local/SINGE/SINGE.sh /usr/local/MATLAB/MATLAB_Runtime/v94 standalone',
                             '/usr/local/SINGE/' + inputMat, '/usr/local/SINGE/' + geneListMat, outFileSymlink, '/usr/local/SINGE/' + paramsFile, '"'])
        print(cmdToRun)
        # also print the parameters
        print("\tParameters: %s" % (', '.join("%s: %s" % (p, str(params[p])) for p in params_order)))
        subprocess.check_call(cmdToRun, shell=True)


def parseOutput(RunnerObj):
    '''
    Function to parse outputs from SINGE.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/SINGE/"
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for idx in range(len(colNames)):  
        
        # Quit if output directory does not exist
        if not Path(outDir+ str(idx)+'/SINGE_Ranked_Edge_List.txt').exists():
            print(outDir+ str(idx)+'/SINGE_Ranked_Edge_List.txt does not exist, skipping...')
            return

        # Read output
        OutSubDF[idx] = pd.read_csv(outDir+ str(idx)+'/SINGE_Ranked_Edge_List.txt',
                            sep = '\t', header = 0)
    # megre the dataframe by taking the maximum value from each DF
    # Code from here: 
    # https://stackoverflow.com/questions/20383647/pandas-selecting-by-label-sometimes-return-series-sometimes-returns-dataframe
    outDF = pd.concat(OutSubDF)
    outDF.columns= ['Gene1','Gene2','EdgeWeight']
    # Group by rows code is from here:
    # https://stackoverflow.com/questions/53114609/pandas-how-to-remove-duplicate-rows-but-keep-all-rows-with-max-value
    res = outDF[outDF['EdgeWeight'] == outDF.groupby(['Gene1','Gene2'])['EdgeWeight'].transform('max')]
    # Sort values in the dataframe   
    finalDF = res.sort_values('EdgeWeight', ascending=False)   
    finalDF.to_csv(outDir+'rankedEdges.csv',sep='\t', index = False)
