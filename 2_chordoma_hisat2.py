#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2019 Charles Lin

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


#scripts for running chordoma slam bams through hisat2




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



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#slam_data_file
slam_data_file = '%sdata_tables/CHORDOMA_CH22_SLAM_TABLE_1.txt' % (projectFolder)

slam_data_file_2 = '%sdata_tables/CHORDOMA_CH22_SLAM_TABLE_2.txt' % (projectFolder)




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
    pipeline_dfci.summary(slam_data_file)

    pipeline_dfci.summary(slam_data_file_2)




    print('\n\n')
    print('#======================================================================')
    print('#===========================II. CALL BOWTIE2===========================')
    print('#======================================================================')
    print('\n\n')

    pipeline_dfci.makeBowtieBashJobsSlurm(slam_data_file_2,namesList = [],launch=True,overwrite=False,pCount=32,paramString=' ')

    # print('\n\n')
    # print('#======================================================================')
    # print('#============================II. WRAP HISAT2===========================')
    # print('#======================================================================')
    # print('\n\n')

    # #quick wrapper for hisat2
    
    # def wrapHisat2(data_file,genome,output_folder,param_string='',launch=True):

    #     '''
    #     simple hisat2 wrapper draws index from a dictionary
    #     '''

    #     genome_dict = {'hg19': '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index/',
    #                    'hg19_ercc':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Hisat2Index_ERCC/'

    #     )

    #     analysis_name = data_file.split('/')[-1].split('.')[0]
    #     hisat_bash_path = '%s%s_hisat2.sh' % (output_folder,analysis_name)

    #     hisat_bash = open(hisat_bash_path,'w')

    #     hisat_bash.write('#!/usr/bin/bash\n\n')
    #     hisat_bash.write('#Running hisat2 on %sn\n' % (analysis_name))

    #     hisat2_cmd = '%s -x %s -U %s -S %s' % (hisat2_path,index_path,

        





#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
