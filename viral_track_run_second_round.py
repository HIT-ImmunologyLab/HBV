#!/usr/bin
#-*-coding:utf-8-*-

import os
import threading
import multiprocessing

def mkdir(dirPath):
    cmd_mkdir = 'mkdir -p %s'%dirPath
    os.system(cmd_mkdir)

def get_virusite_genome_dict(virusite_genome_file,virusite_genome_dict):
    with open(virusite_genome_file,'r')as fin:
        content = fin.read()
        elems = content.split('\n>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            cur_genome_id = elem.split('\n')[0].split('|')[1]
            if cur_genome_id not in virusite_genome_dict:
                if elem.startswith('>'):
                    virusite_genome_dict[cur_genome_id] = elem
                else:
                    virusite_genome_dict[cur_genome_id] = '>'+elem

def get_virus_genome_and_annotation_info(all_virus_genome_file,genome_id):
    with open(all_virus_genome_file,'r')as fin:
        contents = fin.read()
        elems = contents.split('\n>')
        if '' in elems:
            elems.remove('')
        for elem in elems:
            if genome_id in elem:
                if elem.startswith('>'):
                    cur_fasta_info = elem
                else:
                    cur_fasta_info = '>'+elem
                cur_def_info = elem.split('\n')[0]
                cur_genome_id = cur_def_info.split(' ')[0].replace('>','')
                cur_genome_def = ' '.join(cur_def_info.split(' ')[1:])
                cur_virus_name = cur_genome_def.split(',')[0]
                cur_genome_len = str(len(''.join(elem.split('\n')[1:])))
                break
    return [cur_fasta_info,cur_genome_id,cur_genome_len,cur_virus_name,cur_genome_def]

def create_fasta_and_annotation_file(fasta_file,annotation_file,virusite_flag,att_genome_info,exclude_genome_id,virusite_genome_dict,virusite_annotation_file):
    fout = open(fasta_file,'w')
    att_genome_fasta_info = att_genome_info[0]
    fout.write(att_genome_fasta_info+'\n')
    if virusite_flag == 'Yes':
        for key in virusite_genome_dict:
            if key in exclude_genome_id or exclude_genome_id in key:
                continue
            fout.write(virusite_genome_dict[key]+'\n')
    fout.close()

    fout = open(annotation_file,'w')
    header_list = ['Name_sequence','Genome_length','Virus_name','Complete_segment_name']
    fout.write('\t'.join(header_list)+'\n')
    att_genome_annotation_info = att_genome_info[1:]
    fout.write('\t'.join(att_genome_annotation_info)+'\n')
    if virusite_flag == 'Yes':
        with open(virusite_annotation_file,'r')as fin:
            lines = fin.readlines()
            for line in lines[1:]:
                if exclude_genome_id in line:
                    continue
                fout.write(line)
    fout.close()

def create_virus_info_for_create_STAR_index(genome_id,exclude_genome_id,all_virus_genome_file,virusite_genome_file,virusite_annotation_file,genome_save_root_dir):
    virusite_genome_dict = {}
    get_virusite_genome_dict(virusite_genome_file,virusite_genome_dict)

    att_genome_info = get_virus_genome_and_annotation_info(all_virus_genome_file, genome_id)

    genome_with_virusite_index_dir = '%s/%s_and_virusite'%(genome_save_root_dir,genome_id)
    mkdir(genome_with_virusite_index_dir)
    genome_with_virusite_fasta_file = '%s/%s_and_virusite.fasta'%(genome_with_virusite_index_dir,genome_id)
    genome_with_virusite_annotation_file = '%s/%s_and_virusite_annotation.txt'%(genome_with_virusite_index_dir,genome_id)
    virusite_flag = "Yes"
    create_fasta_and_annotation_file(genome_with_virusite_fasta_file, genome_with_virusite_annotation_file, virusite_flag,
                                     att_genome_info, exclude_genome_id,
                                     virusite_genome_dict, virusite_annotation_file)

    genome_with_host_index_dir = '%s/%s_and_host'%(genome_save_root_dir,genome_id)
    mkdir(genome_with_host_index_dir)
    genome_with_host_index_fasta_file = '%s/%s_and_host.fasta' % (genome_with_host_index_dir, genome_id)
    genome_with_host_index_annotation_file = '%s/%s_and_host_annotation.txt' % (genome_with_host_index_dir, genome_id)
    virusite_flag = "No"
    create_fasta_and_annotation_file(genome_with_host_index_fasta_file, genome_with_host_index_annotation_file,
                                     virusite_flag,att_genome_info, exclude_genome_id,
                                     virusite_genome_dict, virusite_annotation_file)

def create_STAR_index_with_host(index_prefix,virus_genome_file,host_genome_file):
    cmd_star = 'STAR --runThreadN 20 --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s  %s'%(index_prefix,virus_genome_file,host_genome_file)
    os.system(cmd_star)

def create_STAR_index_without_host(index_prefix,virus_genome_file):
    cmd_star = 'STAR --runThreadN 20 --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s'%(index_prefix,virus_genome_file)
    os.system(cmd_star)

def create_STAR_index(genome_save_root_dir,star_index_root_dir,sample_name,genome_id):
    index_prefix = '%s/%s/%s_and_host'%(star_index_root_dir,sample_name,genome_id)
    mkdir(index_prefix)
    virus_genome_file = '%s/%s_and_host/%s_and_host.fasta' % (genome_save_root_dir, genome_id, genome_id)
    host_genome_file = '/home/data/viral_track_host_genome/genome/Homo_sapiens.GRCh38.dna.chromosome*.fa'
    # create_STAR_index_with_host(index_prefix, virus_genome_file, host_genome_file)
    muthread = threading.Thread(target=create_STAR_index_with_host, args=(index_prefix, virus_genome_file, host_genome_file,))
    muthread.start()

    indexa_with_virusite_prfix = '%s/%s/%s_and_virusite'%(star_index_root_dir,sample_name,genome_id)
    mkdir(indexa_with_virusite_prfix)
    virus_genome_with_virusite_fasta_file = '%s/%s_and_virusite/%s_and_virusite.fasta' % (genome_save_root_dir, genome_id, genome_id)
    create_STAR_index_without_host(indexa_with_virusite_prfix, virus_genome_with_virusite_fasta_file)

    muthread.join()

def create_parameter_file(parameter_file,result_dir,name_run,index_genome="/home/data/viral_track_database/annotation_index_small",viral_annotation_file="/home/software/Viral-Track/Virusite_annotation_file.txt",thread_num=30,minimal_read_mapped=50):
    fout = open(parameter_file,'w')
    strwrite = 'N_thread=%s\nOutput_directory="%s"\nIndex_genome="%s"\nViral_annotation_file="%s"\nName_run="%s"\nLoad_STAR_module=FALSE\nLoad_samtools_module=FALSE\nLoad_stringtie_module=FALSE\nMinimal_read_mapped=%s\n'\
               %(thread_num,result_dir,index_genome,viral_annotation_file,name_run,minimal_read_mapped)
    fout.write(strwrite)
    fout.close()

def create_target_file(target_file,extract_result):
    fout = open(target_file,'w')
    extract_file_list = os.listdir(extract_result)
    for extract_file_item in extract_file_list:
        cur_extract_file_path = '%s/%s'%(extract_result,extract_file_item)
        fout.write(cur_extract_file_path+'\n')
    fout.close()

def run_viral_tack(parameter_file,target_file):
    command = "Rscript Viral_Track_scanning.R %s %s" % (parameter_file, target_file)
    os.system(command)
    command = "Rscript Viral_Track_transcript_assembly.R %s %s" % (parameter_file, target_file)
    os.system(command)
    command = "Rscript Viral_Track_cell_demultiplexing.R %s %s" % (parameter_file, target_file)
    os.system(command)

def run_viral_track(viral_track_result,seq_save_dir,sample_name,result_dir,sample_merge_genome_prefix,index_genome="viral_track_database/annotation_index_small",viral_annotation_file="Virusite_annotation_file.txt",thread_num=30,minimal_read_mapped=50):
    extract_result = seq_save_dir

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


if __name__=='__main__':
    genome_id_list = ['genome_1', 'genome_2','genome_3']
    sample_name_list = ['sample_1', 'sample_2','sample_3']
    seq_save_dir_list = ['genome_1_path','genome_2_path','genome_3_path']
    viral_track_result = "viral_track_result"
    exclude_genome_id = 'exclude_genome_id'

    for i in range(0,len(genome_id_list)):
        genome_id = genome_id_list[i]
        sample_name = sample_name_list[i]
        all_virus_genome_file = 'all_HBV_B_C_virus_strain.fa'
        virusite_genome_file = 'genomes.fasta'
        virusite_annotation_file = 'Virusite_annotation_file.txt'
        genome_save_root_dir = '%s/STAR_genome_dir'%(viral_track_result)
        seq_save_dir = seq_save_dir_list[i]

        create_virus_info_for_create_STAR_index(genome_id, exclude_genome_id, all_virus_genome_file, virusite_genome_file,
        virusite_annotation_file, genome_save_root_dir)

        star_index_root_dir = "%s/STAR_index"%(viral_track_result)
        create_STAR_index(genome_save_root_dir, star_index_root_dir, sample_name, genome_id)

        run_list =['and_host', 'and_other_virus']
        pool = multiprocessing.Pool(processes=2)
        for run_item in run_list:
            result_root_dir = '%s/result'%(viral_track_result)
            mkdir(result_root_dir)
            result_dir = '%s/%s_with_one_HBV_strain/%s_%s' % (result_root_dir, sample_name, genome_id, run_item)
            mkdir(result_dir)
            sample_merge_genome_prefix = sample_name
            index_genome = '%s/%s/%s_%s' % (star_index_root_dir, sample_name, genome_id, run_item)
            viral_annotation_file = '%s/%s_%s/%s_%s_annotation.txt' % (
            genome_save_root_dir, genome_id, run_item, genome_id, run_item)
            thread_num = 30
            minimal_read_mapped = 50
            pool.apply_async(run_viral_track, (
            viral_track_result,sample_name,result_dir,seq_save_dir, sample_merge_genome_prefix, index_genome, viral_annotation_file, thread_num, minimal_read_mapped,))
        pool.close()
        pool.join()
