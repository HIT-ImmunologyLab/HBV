#!/usr/bin
#-*-coding:utf-8-*-

import os
import multiprocessing

def viral_track_first_round_pipline(sample_name,sample_save_dir,result_root_dir,viral_track_index_dir,viral_track_anno_dir,script_save_dir):
    script = "%s/viral_track_run_first_round.py"%(script_save_dir)
    cmd = "python %s --sample_name %s --sample_save_dir %s --result_root_dir %s --viral_track_index_dir %s --viral_track_anno_dir %s"% \
          (script,sample_name,sample_save_dir,result_root_dir,viral_track_index_dir,viral_track_anno_dir)
    print(cmd)
    os.system(cmd)

def create_parameter_file(sample_save_root_dir, sample_id_list, result_root_dir, viral_track_index_result_root_dir,
                          viral_annotation_save_root_dir, script_save_dir,parameter_file):
    fout = open(parameter_file,'w')
    for sample_id_item in sample_id_list:
        result_dir = "%s/%s"%(result_root_dir,sample_id_item)
        sample_save_dir = "%s/%s"%(sample_save_root_dir,sample_id_item)
        strwrite_list = [sample_id_item,sample_save_dir, result_dir,viral_track_index_result_root_dir,
        viral_annotation_save_root_dir,script_save_dir]
        fout.write('\t'.join(strwrite_list)+'\n')
    fout.close()

def get_parameter_list(parameter_file,parameter_list):
    with open(parameter_file,'r')as fin:
        lines = fin.readlines()
        for line in lines:
            content = line.strip('\n').split('\t')
            parameter_list.append(content)

if __name__=='__main__':
    sample_save_root_dir = 'sequence_sample_save_dir'
    sample_id_list = ['sample_id_1','sample_id_2','sample_id_3','sample_id_4']

    result_root_dir = 'viral_track_result'
    ## place annotation_index_host_and_all_B_C_virus_cdhitest_id_0.95_cov_1.0 to viral_track_database
    viral_track_index_result_root_dir = 'viral_track_database'
    ## place all_HBV_B_C_virus_strain_cdhitest_id_0.95_cov_1.0.txt to HBV_genome
    viral_annotation_save_root_dir = 'HBV_genome'
    parameter_file = '%s/viral_track_first_run_parameter.txt' % result_root_dir

    ## place viral_track_run_first_round.py to script
    script_save_dir = 'script'

    ## step1: create parameter file
    create_parameter_file(sample_save_root_dir, sample_id_list, result_root_dir, viral_track_index_result_root_dir,
                          viral_annotation_save_root_dir, script_save_dir, parameter_file)

    ## step2: get parameter
    parameter_list = []
    get_parameter_list(parameter_file, parameter_list)

    ## step3: run viral track
    ## The maximum number of targets should be assigned based on server performance
    max_task_num = 4
    pool = multiprocessing.Pool(processes=max_task_num)
    for i in range(0, len(parameter_list)):
        task_id = i
        run_parameter = parameter_list[task_id]
        sample_name = run_parameter[0]
        sample_save_dir = run_parameter[1]
        result_root_dir = run_parameter[2]
        viral_track_index_dir = run_parameter[3]
        viral_track_anno_dir = run_parameter[4]
        script_save_dir = run_parameter[5]
        pool.apply_async(viral_track_first_round_pipline, (
            sample_name, sample_save_dir, result_root_dir, viral_track_index_dir, viral_track_anno_dir,
            script_save_dir,))

    pool.close()
    pool.join()



