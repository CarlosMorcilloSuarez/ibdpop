#!/usr/bin/env python3 	
# -*- coding: utf-8 -*-


''' ibdpop

    A tool to work with IBD block information inferred from individuals 
    belonging to different populations
'''

__author__  = "Carlos Morcillo-Suarez"
__license__ = "GPL"
__version__ = "2019/10/04 11:28" # YYYY/MM/DD HH:MM
__email__   = "carlos.morcillo.upf.edu@gmail.com"


import sys
import getopt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def show_df(name,df):
    print()
    print('%s (n=%d)' % (name,len(df)))
    display(df.head())
    print('...')


def lengthInfoExists():
    if 'Length' not in IBDs.columns:
        print()
        print( "** No Length information present in %s" 
                % (ibdsFileName))
        print("** Consider using --impute-length option ")
        sys.exit(2)
    else:
        return(True)
    
    
def processArguments(argv):
    if len(argv) == 0:
        print()
        print("No options given")
        usage()
        sys.exit(2)
    try:
        opts, args = getopt.getopt(
                        argv,
                        "",
                        ["file=","out=","color-code=", 'min-length=',
                        'impute-length',"min-score=","dpi="]
        )
    except getopt.GetoptError as e:
        print(e)
        print('''
        For details of use:
            ibdpop --help 
        ''')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("--file"):
            global dataPrefix
            dataPrefix = arg
        if opt in ("--out"):
            global outputPrefix
            outputPrefix = arg
        if opt in ("--color-code"):
            global colorCode
            colorCode = arg
        if opt in ("--min-length"):
            global minLength
            minLength = int(arg)
        if opt in ("--impute-length"):
            global imputeLength
            imputeLength = True
        if opt in ("--min-score"):
            global minScore
            minScore = int(arg)       
        if opt in ("--dpi"):
            global dpi
            dpi = int(arg)
                        
def usage():
    print('''
        ibdpop
                a tool to analyse IBD segments derived from individuals
                belonging to different populations.
        
        SYNOPSIS
        
            ibdpop <command> [options]
            
        COMMANDS
        
            version            Shows the current version of the program 
            
            help               Shows this help page
        
            matrixcounts       Matrix of the counts of IBDs between all pairs
                               of individuals
                               
            quality            Plots scores vs length of IBDs for QC
            
            recode             Creates a new dataset applying the current 
                               filters
            
        
        OPTIONS
                    
            --file <prefix>
                    Prefix of the files with IBDs (prefix.ibd) and individuals 
                    (prefix.ind) information
                    
            --out <prefix>
                    Prefix that will be used in all files created during the
                    execution
                    
            --color-code <file>
                    File containing the color codes for each population to be
                    used in plots.
                    
            --min-length <length>
                    Filters out IBD segments shorter of <length> cM
                    
            --impute-length
                    Imputes IBD segments length from chromosomal position
                    instead of using provided cM length
                    
            --min-score <score>
                    Filters out IBD segments with score smaller that <score>
                    
            --dpi <dpi>
                    dpi (dots per inch) to be used plots
                    Default = 100
            
    ''')


def plotLengthScore(IBDs,outputPrefix,dpi):
    '''
       Creates plots to evaluate the quality of the list of IBDs
    
       Arguments:       
       IBDs -- pandas dataframe with all pair of IBD segments to plot
       outputPrefix -- prefix of the plot output file       
    '''
    
    lenghtHistogramFileName = outputPrefix+'.lhist.png'
    scoreHistogramFileName = outputPrefix+'.shist.png'
    lengthsVSscoresFileName = outputPrefix+'.lvs.png'
            
    # Length histogram        
    figure = plt.figure(figsize=(12,8),dpi=dpi)
    figure.suptitle("IBDs Length Histogram",fontsize=18)
    bins=int(max(IBDs.Length))*2
    color="blue"

    ax = figure.add_subplot(111)
    ax.hist(IBDs.Length,
            bins=bins,
            edgecolor="black")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    figure.savefig(lenghtHistogramFileName)
    print('Created file: %s' % (lenghtHistogramFileName))
    plt.close()

    # Print scores histogram
    figure = plt.figure(figsize=(12,8),dpi=dpi)
    figure.suptitle("IBDs Score Histogram",fontsize=18)
    bins=int(max(IBDs.Score))
    color="coral"

    ax = figure.add_subplot(111)
    ax.hist(IBDs.Score,
            bins=bins,
            edgecolor="black",
            color=color)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    figure.savefig(scoreHistogramFileName)
    print('Created file: %s' % (scoreHistogramFileName))
    plt.close()

    # Lengths vs scores plot
    figure = plt.figure(figsize=(12,8),dpi=dpi)
    figure.suptitle("IBDs",fontsize=18)
    ax = figure.add_subplot(111)

    ax.plot(IBDs.Length,IBDs.Score,"o",
                markersize = 2)
    ax.set_xlabel("Length",fontsize=16)
    ax.set_ylabel("Score",fontsize=16)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    figure.savefig(lengthsVSscoresFileName)
    print('Created file: %s' % (lengthsVSscoresFileName))    
    plt.close()
    
    
def plotIBDMatrix(IBDs,individuals,outputPrefix,
                  colorOfPopulation = None,dpi=300):
    '''
       Plots the matrix of IBD numbers shared by all pairs of individuals
    
       Arguments:       
       IBDs -- pandas dataframe with all pair of IBD segments to plot
       individuals -- pandas dataframe with population and subpopulation info         
       outputPrefix -- prefix of the plot output file       
       
       Keyword arguments:             
       colorOfPopulation -- dictionary with the '#aabbcc' color code to use for each
                            population 
       dpi -- dpi to be used in the generated plot
    '''

    # Colors of populations
    if not colorOfPopulation:
        colorOfPopulation = {}
        colors=['#5373e8','#9a6324','#f0f69f','#c0c0f0','#4363d8',
                '#ffe119','#b13ed4','#e6beff','#909090','#eeeab8',
                '#808000','#3cb44b','#f59241','#f032e6','#fabebe',
                '#bbffd3','#46f0f0','#911eb4','#e6194b','#ffd8b1',
                '#008080','#808080','#bcf60c','#46f0f0','#a05050',
                '#ffe119','#f58231','#ff497b','#ffd8b1']
        for index,population in enumerate(set(individuals.Pop)):
            colorOfPopulation[population] = colors[index%len(colors)]

    # Creates ordered list of all individuals
    orderedIndividuals = individuals.sort_values(
                                ['Pop','SubPop'],
                                ascending = True)
    
    # Creates matrix
    matrix = pd.DataFrame(index = orderedIndividuals.index, 
                          columns = orderedIndividuals.index)
    matrix.iloc[:,:] = 0

    # Counts IBDs and fill matrix
    for index, IBD in IBDs.iterrows():
        matrix.at[IBD.Ind1,IBD.Ind2] += 1
        matrix.at[IBD.Ind2,IBD.Ind1] += 1
   
    # Defines colorMap
    cmap = plt.cm.winter_r
    cmap.set_under(color='white')    
    
    # Plots Matrix
    height = len(individuals)/8
    width = height * 1.5
    
    fig = plt.figure(figsize = (width,height),
                     dpi=dpi)
    ax = fig.add_subplot(1,1,1)
    img = ax.imshow(matrix,
              cmap=cmap,
              vmin=0.0000001)

    # -- ticks
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.xaxis.set_ticks_position('top')

    ax.set_xticks(range(0,len(matrix.columns)))
    ax.set_xticklabels(matrix.columns,
                       fontsize = 4,
                       rotation = 40 ,
                       horizontalalignment='left',
                       rotation_mode = 'anchor')

    ax.set_yticks(range(0,len(matrix.columns)))
    ax.set_yticklabels(matrix.columns,
                       fontsize = 5)

    # Gridline based on minor ticks
    ax.set_xticks(np.arange(-.5, len(matrix.columns), 1), minor=True);
    ax.set_yticks(np.arange(-.5, len(matrix.columns), 1), minor=True);
    ax.grid(which='minor', color='lightgray', linestyle='-', linewidth=0.1)

    # Grid of subpopulation and subpopulation name
    start = -1
    for index in range(len(matrix.columns)-1):
        if individuals.at[matrix.columns[index],'SubPop'] \
           != individuals.at[matrix.columns[index+1],'SubPop']:
            ax.axhline(index+0.5,
                      linewidth = 0.2)
            ax.axvline(index+0.5,
                      linewidth = 0.2)
            ax.text(-8,
                    start+1+(index-start)/2,
                    individuals.at[matrix.columns[index],'SubPop'],
                    fontsize = 6,
                    weight = 'bold',
                    horizontalalignment = 'right')
            start = index
            
    ax.text(-8,
            start+1+(len(matrix.columns)-start)/2,
            individuals.at[matrix.columns[-1],'SubPop'],
            fontsize = 6,
            weight = 'bold',
            horizontalalignment = 'right')

    # Grid of population and population name
    start = 0
    for index in range(len(matrix.columns)-1):
        if individuals.at[matrix.columns[index],'Pop'] \
           != individuals.at[matrix.columns[index+1],'Pop']:
            ax.axhline(index+0.5,
                      linewidth = 1,
                      color = 'black')
            ax.axvline(index+0.5,
                      linewidth = 1,
                      color = 'black')
            ax.text(len(matrix.columns)+2,
                    start+(index-start)/2,
                    individuals.at[matrix.columns[index],'Pop'])
            start = index + 1
    ax.text(len(matrix.columns)+2,
            start+(len(matrix.columns)-start)/2,
            individuals.at[matrix.columns[-1],'Pop'])

    # Colors and names of marks of population print('Created file: %s' % (outputPrefix+'.cmatrix.png'))
    for index in range(len(matrix.columns)):
        color = colorOfPopulation[individuals.at[matrix.columns[index],'Pop']]
        ax.add_patch(
            Rectangle((len(matrix.columns)-0.5,index-0.5),2,1,                 
                       clip_on=False,
                       color = color))

    # Adds Colorbar
    ax_colorbar = fig.add_axes([0,0,0.05,1]) 
    ax_colorbar.axis("off")
    fig.colorbar(img,ax=ax_colorbar,shrink = 0.2)     

    fig.savefig(outputPrefix+'.cmatrix.png')
    
    print('Created file: %s' % (outputPrefix+'.cmatrix.png'))


def loadIBDs(ibdsFileName,imputeLength = False):
    '''
    Returns a df with the contents of the ibdsFileName.
    If imputeLength == True, Length is calculated from chromosome positions
    '''
    IBDs = pd.read_csv(ibdsFileName,sep='\t')
    if imputeLength:
        IBDs["Length"] = (IBDs['End'] - IBDs['Start']) / 1000000
    return(IBDs)


if __name__ == "__main__":
    
    print('ibdpop - Version: %s' % (__version__))
    print('(C) 2019 Carlos Morcillo-Suarez   GPL License')
    print()

    # Reads command
    if len(sys.argv) < 2 or sys.argv[1] not in ['version',
                                                'help',
                                                'matrixcounts',
                                                'quality',
                                                'recode']:
        print('No command provided - EXITING')
        usage()
        sys.exit(2)
    else:
        command = sys.argv[1]

    # Global Variables
    dataPrefix = ""
    outputPrefix = "ibdpop"
    colorCode = ""
    minLength = 0
    minScore = 0
    imputeLength = False
    dpi = 100

    # version
    if command == 'version':
        print('ibdpop version: '+__version__)
        sys.exit()
        
    # help
    if command == 'help':
        usage()
        sys.exit() 

    # Process command line
    processArguments(sys.argv[2:])
    
    # Uploads individuals
    individualsFileName = dataPrefix+".ind"
    individuals = pd.read_csv(individualsFileName,sep='\t',index_col=0)
    print("Uploaded %d individuals from %s" 
            % (len(individuals),
               individualsFileName) )

    # Uploads IBDs
    ibdsFileName = dataPrefix+".ibd"
    IBDs = loadIBDs(ibdsFileName,imputeLength=imputeLength)
    print("Uploaded %d IBD segments from %s" 
            % (len(IBDs),
               ibdsFileName) )
    if minLength > 0:
        if lengthInfoExists():
            IBDs = IBDs[IBDs.Length > minLength]
    if minScore > 0:
        IBDs = IBDs[IBDs.Score > minScore]
    print("%d final IBD segments after length and score filters" 
          % (len(IBDs)))    
        
    # Uploads color code of populations
    colorOfPopulation = {}
    if colorCode:
        print("Reads color codes from %s" % (colorCode))
        for index, row in pd.read_csv(colorCode,sep='\t').iterrows():
            colorOfPopulation[row.Population] = row.Color
   
   
    # COMMANDS ---------------------------------------
    # matrixcounts       
    if command == 'matrixcounts': 
        plotIBDMatrix(IBDs,
                      individuals,
                      outputPrefix,
                      colorOfPopulation = colorOfPopulation,
                      dpi=dpi)
                      
    # quality          
    elif command == 'quality':
        if lengthInfoExists():
            plotLengthScore(IBDs, outputPrefix,dpi)
        
    # recode
    elif command == 'recode':
        IBDs.to_csv(outputPrefix+'.ibd', sep='\t', index = False)
        print('Created file: %s' % (outputPrefix+'.ibd'))    
        individuals.to_csv(outputPrefix+'.ind', sep='\t')
        print('Created file: %s' % (outputPrefix+'.ind'))    
    
    
    
    
