#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of slam seq analysis from Muhar et al., 2018




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'chordoma_slam'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files


#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)


py27_path = '/storage/cylin/anaconda3/envs/py27_anaconda/bin/python'
#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
chip_data_file = '%sdata_tables/CHORDOMA_CH22_CHIP.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#===========================II. CALLING MACS===========================')
    print('#======================================================================')
    print('\n\n')

    #pipeline_dfci.run_macs(chip_data_file,projectFolder,macsFolder,macsEnrichedFolder,wiggleFolder,useBackground=True)


    

    print('\n\n')
    print('#======================================================================')
    print('#===================II. RUNNING ENHANCER PROMOTER ON T=================')
    print('#======================================================================')
    print('\n\n')


    def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list = [],useBackground=True):

        '''
        runs enhancer promoter on everybody with the conserved regions and union of active genes
        '''

        #hard coded paths
        tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

        #setting the output folder
        ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)

        dataDict = pipeline_dfci.loadDataTable(dataFile)
        if len(names_list) == 0:
            names_list = [name for name in dataDict.keys()]
            names_list.sort()

        bams_list = [dataDict[name]['bam'] for name in names_list]
        bams_string = ' '.join(bams_list)

        background_names = [dataDict[name]['background'] for name in names_list]
        background_list = [dataDict[background_name]['bam'] for background_name in background_names]
        background_string = ' '.join(background_list)


        ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
        ep_bash = open(ep_bash_path,'w')

        ep_bash.write('#!/usr/bin/bash\n\n\n')

        ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))

        if useBackground:
            python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

            ep_bash.write(python_cmd)

        else:
            python_cmd = 'python %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

            ep_bash.write(python_cmd)

        ep_bash.close()

        return(ep_bash_path)



    active_gene_path = '%sgeneListFolder/HG19_CH22_ACTIVE.txt' % (projectFolder)
    t_path = '%sCH22_dTag_T_WT_HA.bed' % (macsEnrichedFolder)
    analysis_name = 'CH22_T'
    bam_list = ['CH22_dTag_T_WT_HA']    
    wrap_enhancer_promoter(chip_data_file,t_path,active_gene_path,analysis_name,bam_list,useBackground=True)

        

    print('\n\n')
    print('#======================================================================')
    print('#=======================IV. PLOTTING REGIONS===========================')
    print('#======================================================================')
    print('\n\n')

    figure_gff = [['chr7','CAV1','CAV1',116131753,116236016,'','+','','CAV1'],
                  ['chr12','KRT18','KRT18',53286891,53351373,'','+','','KRT18'],
                  ['chr1','MCL1','MCL1',150530141,150556155,'','+','','MCL1'],
                  ]


    figure_gff_path = '%sHG19_CHORDOMA_SLAM_FIGURES.gff' % (gffFolder)
    utils.unParseTable(figure_gff,figure_gff_path,'\t')

    plotName = 'HG19_CHORDOMA_SLAM_FIGURES'
    outputFolder = utils.formatFolder('%sgene_plot' % (projectFolder),True)
    plot_list = ['CH22_dTag_T_WT_HA','CH22_H3K27AC']
    pipeline_dfci.callBatchPlot(chip_data_file,figure_gff_path,plotName,outputFolder,plot_list,uniform=False,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '',scaleFactorString ='')

                  




#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
