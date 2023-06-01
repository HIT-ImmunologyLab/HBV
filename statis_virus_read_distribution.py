#!/usr/bin
#-*-coding:utf-8-*-

import os
import argparse

def statis_mapped_coverage(bam_file_path,cov_file_path):
    cmd_bedtools = 'bedtools genomecov -dz -ibam "%s" > "%s"'%(bam_file_path,cov_file_path)
    os.system(cmd_bedtools)

def plot_virus_coverage(script_save_dir,cov_file,result_file):
    cmd_script = "/usr/bin/Rscript %s/plot_virus_read_distribution.R %s %s"%(script_save_dir,cov_file,result_file)
    os.system(cmd_script)

def virus_mapping_analysis(bam_file,cov_file,virus_plot_file,script_save_dir):
    statis_mapped_coverage(bam_file, cov_file)
    plot_virus_coverage(script_save_dir, cov_file, virus_plot_file)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_file', help="bam file of mapping virus")
    parser.add_argument('--cov_file', help="coverage file of mapping virus")
    parser.add_argument('--virus_plot_file', help="plot file of virus mapping coverage")
    parser.add_argument('--script_dir', help="directory of script")

    args = parser.parse_args()
    if args.bam_file:
        bam_file = args.bam_file
    if args.cov_file:
        cov_file = args.cov_file
    if args.virus_plot_file:
        virus_plot_file = args.virus_plot_file
        if args.script_dir:
            script_dir = args.script_dir

    virus_mapping_analysis(bam_file, cov_file, virus_plot_file,script_dir)
