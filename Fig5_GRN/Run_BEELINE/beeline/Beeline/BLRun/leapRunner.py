import os
import pandas as pd
from pathlib import Path
from .singularityUtils import get_singularity_image_path
import numpy as np

def generateInputs(RunnerObj):
    '''
    If the folder/files under RunnerObj.datadir exist, 
    this function will not do anything.
    '''
    if not RunnerObj.inputDir.joinpath("LEAP").exists():
        print("Input folder for LEAP does not exist, creating input folder...")
        RunnerObj.inputDir.joinpath("LEAP").mkdir(parents=True, exist_ok=True)
        
    ExpressionData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.exprData),
                                     header = 0, index_col = 0)
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        # Select cells belonging to each pseudotime trajectory
        colName = colNames[idx]
        index = PTData[colName].index[PTData[colName].notnull()]
        exprName = "LEAP/ExpressionData"+str(idx)+".csv"
        
        subPT = PTData.loc[index,:]
        subExpr = ExpressionData[index]
        # Order columns by PseudoTime
        newExpressionData = subExpr[subPT.sort_values([colName]).index.astype(str)]
        
        newExpressionData.insert(loc = 0, column = 'GENES', \
                                                     value = newExpressionData.index)


        # Write .csv file
        newExpressionData.to_csv(RunnerObj.inputDir.joinpath(exprName),
                             sep = ',', header  = True, index = False)
    
def run(RunnerObj):
    '''
    Function to run LEAP algorithm

    Requires the maxLag parameter
    '''
    
    inputPath = "data" + str(RunnerObj.inputDir).split(str(Path.cwd()))[1]
    
    maxLag = str(RunnerObj.params['maxLag'])
    
    # make output dirs if they do not exist:
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/LEAP/"
    os.makedirs(outDir, exist_ok = True)
    
    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    for idx in range(len(colNames)):
        exprName = "/LEAP/ExpressionData"+str(idx)+".csv"
        outPath = 'data/' +  str(outDir) + 'outFile'+str(idx)+'.txt'

        # Convert to Singularity with proper bind mount and container execution
        singularity_image = get_singularity_image_path('grnbeeline/leap:base')
        
        cmdToRun = ' '.join(['singularity exec --bind', 
                             str(Path.cwd())+':/data', 
                             str(singularity_image), '/bin/sh -c "mkdir -p /data/' + str(outDir) + 
                             ' && time -v -o /data/' + str(outDir) + 'time'+str(idx)+'.txt', 
                             'Rscript /runLeap.R',
                             '/' + inputPath + exprName, maxLag, '/data/' + str(outDir) + 'outFile'+str(idx)+'.txt', '"'])
        print(cmdToRun)
        os.system(cmdToRun)



def parseOutput(RunnerObj):
    '''
    Function to parse outputs from LEAP.
    '''
    outDir = "outputs/"+str(RunnerObj.inputDir).split("inputs/")[1]+"/LEAP/"

    PTData = pd.read_csv(RunnerObj.inputDir.joinpath(RunnerObj.cellData),
                             header = 0, index_col = 0)

    colNames = PTData.columns
    OutSubDF = [0]*len(colNames)

    for indx in range(len(colNames)):
        outFileName = 'outFile'+str(indx)+'.txt'
        # Quit if output file does not exist
        if not Path(outDir+outFileName).exists():
            print(outDir+outFileName+' does not exist, skipping...')
            return
        
        # Read output
        OutSubDF[indx] = pd.read_csv(outDir+outFileName, sep = '\t', header = 0)
        OutSubDF[indx].Score = np.abs(OutSubDF[indx].Score)
    outDF = pd.concat(OutSubDF)
    FinalDF = outDF[outDF['Score'] == outDF.groupby(['Gene1','Gene2'])['Score'].transform('max')]

    outFile = open(outDir + 'rankedEdges.csv','w')
    outFile.write('Gene1'+'\t'+'Gene2'+'\t'+'EdgeWeight'+'\n')

    for idx, row in FinalDF.sort_values(['Score'], ascending = False).iterrows():
        outFile.write('\t'.join([row['Gene1'],row['Gene2'],str(row['Score'])])+'\n')
    outFile.close()
    
