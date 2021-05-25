65;5604;1c#!/usr/bin/python


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


#scripts for chordoma slam seq analysis




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
    #pipeline_dfci.summary(slam_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#======================II. SET UP SLAM TABLE===========================')
    print('#======================================================================')
    print('\n\n')


    # #for original pilot experiment
    # slam_parent_folder = utils.formatFolder('%sslam/' % (projectFolder),True)
    # slam_table_path = '%sslam_table.txt' % (slam_parent_folder)
    # create_slam_table(slam_data_file,names_list = [],output = slam_table_path)

    # analysis_name = 'ch22_slam'
    # slam_bash_path = '%sch22_slam.sh' % (slam_parent_folder)
    # wrap_slam(genome,analysis_name,slam_table_path,slam_parent_folder,bed_type='utr',slam_bash_out=slam_bash_path)


    # slam_parent_folder = utils.formatFolder('%sslam_190412/' % (projectFolder),True)
    # slam_table_path = '%sslam_table.txt' % (slam_parent_folder)
    # create_slam_table(slam_data_file_2,names_list = [],output = slam_table_path)

    # analysis_name = 'ch22_slam'
    # slam_bash_path = '%sch22_slam.sh' % (slam_parent_folder)
    # wrap_slam(genome,analysis_name,slam_table_path,slam_parent_folder,bed_type='utr',slam_bash_out=slam_bash_path)


    print('\n\n')
    print('#======================================================================')
    print('#======================II. SET UP SLAM TABLE===========================')
    print('#======================================================================')
    print('\n\n')


    slam_parent_folder = utils.formatFolder('%sslam_190412/' % (projectFolder),True)

    slam_table_path = '%sslam_table.tsv' % (slam_parent_folder)
    analysis_name = 'ch22_slam'
    counts_folder = '%s%s_slamdunk/count/' % (slam_parent_folder,analysis_name)
    formatSlamCounts(annotFile,slam_table_path,counts_folder,[],analysis_name,output_folder= '')

#==========================================================================
#==========================ADDITIONAL FUNCTIONS============================
#==========================================================================



def create_slam_table(data_file,names_list= [],output=''):

    '''
    creates a tsv formatted for slam-seq analysis that incorporates time series info
    '''

    data_dict = pipeline_dfci.loadDataTable(data_file)

    if names_list == []:
        names_list = data_dict.keys()
    names_list.sort()

    slam_tsv = []
    for name in names_list:
        fastq_path = data_dict[name]['fastq']
        timepoint = [element for element in name.split('_') if element.count('hr') == 1]
        if len(timepoint) == 0:
            timepoint ='1'
        else:
            timepoint = timepoint[0].replace('hr','')
        new_line = [fastq_path,name,name,timepoint]
        slam_tsv.append(new_line)

    if output != '':

        if utils.checkOutput(output,0,0):
            print('slam table found at %s' % (output))
            print('not overwriting output')
            return output
        else:
            print('writing slam table to %s' % (output))
            utils.unParseTable(slam_tsv,output,'\t')
            return output
    else:
        return slam_tsv



def wrap_slam(genome,analysis_name,slam_table_path,slam_parent_folder,bed_type='utr',slam_bash_out=''):

    '''
    creates a bash file to run slam_seq w/ default parameters on a slam table
    '''
    slam_table = utils.parseTable(slam_table_path,'\t')
    bam_name_list = [line[0].split('/')[-1].replace('.gz','') for line in slam_table]
    print(bam_name_list)

    slam_output_folder = utils.formatFolder('%s%s_slamdunk' % (slam_parent_folder,analysis_name),True)
    separated_reads_folder = utils.formatFolder('%sfilter_split' % (slam_output_folder),True)
    #set the slam bash output file
    if slam_bash_out == '':
        #make up a name based off the analysis_name
        slam_bash_out = '%s%s_slamdunk_cmd.sh' % (slam_parent_folder,analysis_name)

    #get the genome stuff sorted out
    genome_dict = {'hg38':{'genome_fasta':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome_ercc.fa',
                           'whole_gene_bed':'/storage/thinc/projects/dhx15/beds/hg38_RefSeq_Curated_Whole_Gene_ERCC.bed',
                           'utr_bed':''
                           },
                   'hg19':{'genome_fasta':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_ercc.fa',
                           'whole_gene_bed':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/beds/hg19_RefSeq_Curated_Whole_Gene_ercc.bed',
                           'utr_bed':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/beds/hg19_RefSeq_Curated_3UTR_ercc.bed'
                           }
                   }

    genome_fasta = genome_dict[genome]['genome_fasta']

    if bed_type == 'whole':
        bed_path = genome_dict[genome]['whole_gene_bed']
    elif bed_type == 'utr':
        bed_path = genome_dict[genome]['utr_bed']
    else:
        print('ERROR: SPECIFY A BED TYPE {whole,utr}')
        sys.exit()

    slam_cmd = 'slamdunk all -r %s -b %s -o %s -5 12 -n 100 -t 32 -m -rl 75 --skip-sam %s' % (genome_fasta,bed_path,slam_output_folder,slam_table_path)

    alleyoop_cmd = '#alleyoop read-separator -o %s -s %ssnp/ -r %s -t 32 %sfilter/*no_ercc.bam' % (separated_reads_folder,slam_output_folder,genome_fasta,slam_output_folder)

    bash_fh = open(slam_bash_out,'w')
    bash_fh.write('#!/usr/bin/bash\n\n')
    bash_fh.write('source activate slamdunk\n\n')
    #bash_fh.write('#SBATCH --cpus-per-task=32\n')
    #bash_fh.write('#SBATCH --mem=500000\n')
    #bash_fh.write('#SBATCH -n=8\n')
    bash_fh.write(slam_cmd+'\n\n')
    bash_fh.write('cd %sfilter/\n\n'% (slam_output_folder))
    #remove erccs from all of the bams
    for bam_name in bam_name_list:
        bam_path = '%sfilter/%s_slamdunk_mapped_filtered.bam' % (slam_output_folder,bam_name)
        out_bam_path = bam_path.replace('.bam','.no_ercc.bam')
        view_cmd = 'samtools view -h %s | grep -v "ERCC-" | samtools view -Sb - > %s' % (bam_path,out_bam_path)
        index_cmd = 'samtools index %s' % (out_bam_path)
        bash_fh.write(view_cmd + '\n')
        bash_fh.write(index_cmd + '\n')



    bash_fh.write(alleyoop_cmd + '\n\n')
    bash_fh.close()

    print('wrote slamdunk cmd to %s' % (slam_bash_out))
    return slam_bash_out



def formatSlamCounts(annotFile,slam_table_path,counts_folder,names_list =[],analysis_name='',output_folder= ''):

    '''
    summarizes slam dunk output to produce a table of CPM and conversion fraction
    '''

    if analysis_name =='':
        analysis_name = slam_table_path.split('/')[-1].split('.')[0]

    print('Running analysis on %s' % (analysis_name))

    if output_folder == '':
        output_folder= utils.getParentFolder(slam_table_path)
        print('No output folder specified. Writing output to %s' % (output_folder))
    else:
        output_folder = utils.formatFolder(output_folder,True)
        print('writing output to %s' % (output_folder))

    counts_folder = utils.formatFolder(counts_folder,False)

    print('Making start dict from %s' % (annotFile))
    start_dict = utils.makeStartDict(annotFile)
 
    print('Loading slam annotation table from %s:' % (slam_table_path))
    slam_table = utils.parseTable(slam_table_path,'\t')
    
    slam_dict = {}
    for line in slam_table:
        slam_dict[line[1]] = line[0].split('/')[-1].split('.')[0]

    if len(names_list) == 0:
        names_list = slam_dict.keys()
        names_list.sort()

    print(names_list)

    #setting up the header
    header = ['REF_ID','GENE_NAME']
    for name in names_list:
        header += ['%s_CONVERSION' % (name),'%s_CPM' % (name)]

    print(header)

    print('Finding gene names for count table')
    raw_counts_table = [header]
    #use the first one to get the name/refid
    count_table = utils.parseTable('%s%s.fastq_slamdunk_mapped_filtered_tcount.tsv' % (counts_folder,slam_dict[names_list[0]]),'\t')
    for line in count_table[3:]:
        ref_id = line[3].split('.')[0]
        if ref_id in start_dict:
            name = start_dict[ref_id]['name']
        else:
            name = 'NA'
        raw_counts_table.append([ref_id,name])


    #now we go through each dataset
    for name in names_list:
        print('processing %s' % name)
        count_table = utils.parseTable('%s%s.fastq_slamdunk_mapped_filtered_tcount.tsv' % (counts_folder,slam_dict[name]),'\t')
        for i in range(3,len(count_table)):
            line = count_table[i]
            table_index = i-2
            # print(i)
            # print(line)
            # print(table_index)
            # print(raw_counts_table[table_index])

            raw_counts_table[table_index]+= [round(float(line[6]),5),round(float(line[7]),5)]

    #now we need to make a filtered table where we take the highest avg cpm per gene
    #keyed by gene name
    counts_dict = defaultdict(list)
    avg_dict = defaultdict(list)
    for line in raw_counts_table[1:]:
        if line[1] == 'NA':
            continue
        counts_dict[line[1]].append(line)
        cpm_line = [line[n] for n in range(3,len(line),2)]
        avg_dict[line[1]].append(numpy.mean(cpm_line))

    #now get the gene list
    gene_list = counts_dict.keys()
    gene_list.sort()

    filtered_table = [header]
    for gene in gene_list:
        max_cpm = max(avg_dict[gene])
        if max_cpm == 0:
            continue
        #take the first entry with the max cpm
        max_index_list = [i for i in range(len(avg_dict[gene])) if avg_dict[gene][i] == max_cpm]
        ref_string = ','.join([counts_dict[gene][i][0] for i in max_index_list])
        filtered_line = [ref_string] + counts_dict[gene][max_index_list[0]][1:]
        filtered_table.append(filtered_line)


    #now write out to disk
    raw_out_path ='%s%s_RAW.txt' % (output_folder,analysis_name)
    filtered_out_path ='%s%s_FILTERED.txt' % (output_folder,analysis_name)

    utils.unParseTable(raw_counts_table,raw_out_path,'\t')
    utils.unParseTable(filtered_table,filtered_out_path,'\t')


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()



