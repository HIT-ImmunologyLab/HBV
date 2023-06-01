#!/usr/bin
#-*-coding:utf-8-*-

import os
import threading
import multiprocessing
import argparse

def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%(dirPath)
    os.system(cmd_mkdir)

def demultiplexing(sample_genome_save_dir, sample_prefix, demultiplexing_result,extract_result):
    # Step 2: Identify correct cell barcodes
    list_file = os.path.join(demultiplexing_result, 'whitelist.txt')
    command = "umi_tools whitelist --stdin %s/%s_R1.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --set-cell-number=10000 --log2stderr > %s" % (
    sample_genome_save_dir, sample_prefix, list_file)
    os.system(command)

    # Step 3: Extract barcdoes and UMIs and add to read names
    command = "umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin %s/%s_R1.fastq.gz --stdout %s/%s_R1_extract.fastq --read2-in %s/%s_R2.fastq.gz --read2-out=%s/%s_R2_extract.fastq --filter-cell-barcode --whitelist=%s" % (
    sample_genome_save_dir, sample_prefix, extract_result, sample_prefix, sample_genome_save_dir, sample_prefix, extract_result, sample_prefix, list_file)
    os.system(command)


def run_viral_tack(parameter_file,target_file):
    command = "Rscript Viral_Track_scanning.R %s %s" % (
    parameter_file, target_file)
    os.system(command)
    command = "Rscript Viral_Track_transcript_assembly.R %s %s" % (
    parameter_file, target_file)
    os.system(command)
    command = "Rscript Viral_Track_cell_demultiplexing.R %s %s" % (
    parameter_file, target_file)
    os.system(command)

def create_parameter_file(parameter_file,result_dir,name_run,index_genome="viral_track_database/annotation_index_small",viral_annotation_file="Virusite_annotation_file.txt",thread_num=30,minimal_read_mapped=50):
    fout = open(parameter_file,'w')
    strwrite = 'N_thread=%s\nOutput_directory="%s"\nIndex_genome="%s"\nViral_annotation_file="%s"\nName_run="%s"\nLoad_STAR_module=FALSE\nLoad_samtools_module=FALSE\nLoad_stringtie_module=FALSE\nMinimal_read_mapped=%s\n'\
               %(thread_num,result_dir,index_genome,viral_annotation_file,name_run,minimal_read_mapped)
    fout.write(strwrite)
    fout.close()

def create_target_file(target_file,extract_result):
    fout = open(target_file,'w')
    extract_file_list = os.listdir(extract_result)
    for extract_file_item in extract_file_list:
        # if "R2" in extract_file_item:
        cur_extract_file_path = '%s/%s'%(extract_result,extract_file_item)
        fout.write(cur_extract_file_path+'\n')
    fout.close()

def merge_sample_genomes(sample_genome_save_dir,sample_merge_genome_save_dir,sample_merge_genome_prefix):
    command1 = "cat %s/*_R1* > %s/%s_R1.fastq.gz" % (sample_genome_save_dir, sample_merge_genome_save_dir, sample_merge_genome_prefix)
    command2 = "cat %s/*_R2* > %s/%s_R2.fastq.gz" % (sample_genome_save_dir, sample_merge_genome_save_dir, sample_merge_genome_prefix)
    os.system(command1)
    os.system(command2)

def run_viral_track(sample_genome_save_dir,result_dir,sample_merge_genome_prefix,index_genome="/viral_track_database/annotation_index_small",viral_annotation_file="/home/software/Viral-Track/Virusite_annotation_file.txt",thread_num=30,minimal_read_mapped=50):
    prepare_dir = '%s/prepare_for_viral_track'%('/'.join(result_dir.split("/")[0:-1]))
    # merge sample
    sample_merge_genome_save_dir = '%s/sample_merge_genome'%(prepare_dir)
    mkdir(sample_merge_genome_save_dir)
    merge_sample_genomes(sample_genome_save_dir, sample_merge_genome_save_dir, sample_merge_genome_prefix)

    # demultiplexing sample
    demultiplexing_result = '%s/sample_demultiplexing' % (prepare_dir)
    mkdir(demultiplexing_result)
    extract_result = '%s/sample_extract_genome' % (prepare_dir)
    mkdir(extract_result)
    demultiplexing(sample_merge_genome_save_dir, sample_merge_genome_prefix, demultiplexing_result,extract_result)

    # # create parameter file and target file
    parameter_file = '%s/viral_track_parameter.txt'%(result_dir)
    target_file =  '%s/viral_track_target.txt'%(result_dir)
    name_run = sample_merge_genome_prefix
    muthread1 = threading.Thread(target=create_target_file, args=(target_file, extract_result,))
    muthread1.start()
    create_parameter_file(parameter_file, result_dir, name_run,index_genome,viral_annotation_file,thread_num,minimal_read_mapped)
    muthread1.join()

    # run viral-track
    run_viral_tack(parameter_file, target_file)

def viral_track_first_pipeline(sample_name,sample_genome_save_dir,result_root_dir,viral_track_index_dir,viral_track_anno_dir):
    identity_list = ['0.95']
    for identity_item in identity_list:
        result_dir = '%s/merge_read_file_result_host_and_all_B_C_virus_cdhitest_id_%s_cov_1.0' % (
        result_root_dir, identity_item)
        mkdir(result_dir)
        sample_merge_genome_prefix = sample_name
        index_genome = "%s/annotation_index_host_and_all_B_C_virus_cdhitest_id_%s_cov_1.0" % (
        viral_track_index_dir, identity_item)
        viral_annotation_file = "%s/all_HBV_B_C_virus_strain_cdhitest_id_%s_cov_1.0.txt" % (
        viral_track_anno_dir, identity_item, identity_item)
        thread_num = 30
        minimal_read_mapped = 50
        run_viral_track(sample_genome_save_dir, result_dir, sample_merge_genome_prefix, index_genome,
                        viral_annotation_file, thread_num, minimal_read_mapped)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', help="sample name")
    parser.add_argument('--sample_save_dir', help="root directory of samples fastq file")
    parser.add_argument('--result_root_dir', help="directory of result file")
    parser.add_argument('--viral_track_index_dir', help="directory of viral track index")
    parser.add_argument('--viral_track_anno_dir', help="directory of viral track annotation")

    args = parser.parse_args()
    if args.sample_name:
        sample_name = args.sample_name
    if args.sample_save_dir:
        sample_genome_save_dir = args.sample_save_dir
    if args.result_root_dir:
        result_root_dir = args.result_root_dir
    if args.viral_track_index_dir:
        viral_track_index_dir = args.viral_track_index_dir
    if args.viral_track_anno_dir:
        viral_track_anno_dir = args.viral_track_anno_dir

    viral_track_first_pipeline(sample_name, sample_genome_save_dir, result_root_dir, viral_track_index_dir, viral_track_anno_dir)

