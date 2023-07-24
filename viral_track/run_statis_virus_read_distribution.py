#!/usr/bin
#-*-coding:utf-8-*-

import os
import multiprocessing

def virus_mapping_pipline(bam_file,cov_file,virus_plot_file,plot_script_save_dir,script_save_dir):
    script = "%s/statis_virus_read_distribution.py"%(script_save_dir)
    cmd = "python %s --bam_file %s --cov_file %s --virus_plot_file %s --script_dir %s"% \
          (script,bam_file,cov_file,virus_plot_file,plot_script_save_dir)
    print(cmd)
    os.system(cmd)


def create_parameter_file(result_root_dir, sample_id_list,plot_script_save_dir,script_save_dir,parameter_file):
    fout = open(parameter_file,'w')
    for i in range(0,len(sample_id_list)):
        cur_sample_id = sample_id_list[i]
        cur_bam_file = "%s/%s/%s.bam"% (result_root_dir,cur_sample_id,cur_sample_id)
        cur_cov_file = "%s/%s/%s_mapping_genome_cov.txt"%(result_root_dir,cur_sample_id,cur_sample_id)
        cur_virus_plot_file = "%s/%s/%s_mapping_genome_cov.pdf"%(result_root_dir,cur_sample_id,cur_sample_id)
        strwrite_list = [cur_bam_file,cur_cov_file,cur_virus_plot_file,plot_script_save_dir, script_save_dir]
        fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

def get_parameter_list(parameter_file,parameter_list):
    with open(parameter_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            parameter_list.append(content)


def run_viral_track_mapping(result_root_dir,sample_id_list,script_save_dir,plot_script_save_dir,max_task_num):

    parameter_file = '%s/virus_mapping_coverage_parameter.txt' % result_root_dir

    ## step1: create parameter file
    create_parameter_file(result_root_dir, sample_id_list, plot_script_save_dir, script_save_dir,
                          parameter_file)

    ## step2: get parameter
    parameter_list = []
    get_parameter_list(parameter_file, parameter_list)

    ## step3: run cellranger count

    pool = multiprocessing.Pool(processes=min(max_task_num, len(parameter_list)))
    for i in range(0, len(parameter_list)):
        task_id = i
        run_parameter = parameter_list[task_id]
        bam_file = run_parameter[0]
        cov_file = run_parameter[1]
        virus_plot_file = run_parameter[2]
        plot_script_save_dir = run_parameter[3]
        script_save_dir = run_parameter[4]

        pool.apply_async(virus_mapping_pipline,
                         (bam_file, cov_file, virus_plot_file, plot_script_save_dir, script_save_dir,))
    pool.close()
    pool.join()

if __name__=='__main__':
    result_root_dir = 'viral_track_result'
    sample_id_list = ['sample_id1','sample_id2','sample_id3']

    ## place statis_virus_read_distribution.py to analysis/script
    script_save_dir = 'analysis/script'
    ## place plot_virus_read_distribution.R to plot/script
    plot_script_save_dir = 'plot/script'

    ## The maximum number of targets should be assigned based on server performance
    max_task_num = 4

    run_viral_track_mapping(result_root_dir, sample_id_list, script_save_dir, plot_script_save_dir, max_task_num)

