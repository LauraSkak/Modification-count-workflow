######################################################################################
#A.SET ENVIRONMENT AND IMPORT SYS,OS
###################################################################################### 
from gwf import Workflow
import sys, os

gwf = Workflow()

sys.path.append("MODULES")

import basecall_module as bctp
import alignment_module as altp
import quality_control_module as qctp
import variant_call_module as vctp
import modification_count_module as mtp
import util_module as util


######################################################################################
#B.SET CONSTANTS
######################################################################################

### Directories ###

project_path = None # FIXME should be the 

data_path = f'{project_path}/DATA'
results_path = f'{project_path}/RESULTS'

### The workflow parameters ###

parameters = {
    "dorado_model": None, # FIXME should be full path to dorado model
    "remora_model": None, # FIXME should be full path to remora model
    "clair3_model": None, # FIXME if you want phasing, should be full path to clair3 model
    "reference": None, # FIXME should be full path to reference genome
}

mod_name = None # FIXME. This should be the modification name you get on the outfiles.
ref_name = None # FIXME. This should be the reference genome name you get on the outfiles.

sample_dict = util.create_sample_dict(f'{data_path}/sample_meta.tsv') # FIXME. This file should exist and should have a row for each sample with a SOMAseq ID and the desired sample name.

######################################################################################
#C.WORKFLOW
######################################################################################

summary_file_list_pycoQC = []
alignment_file_list_pycoQC = []

gwf.target_from_template(
    name = f"Index_{ref_name}_reference_to_fai", 
    template = altp.fai_index_reference(
        reference = parameters[f"reference"], 
        fai_reference = f'{parameters[f"reference"]}.fai'))

for sample_list in sample_dict.values():
    
    sample_name = sample_list[0]

    basecall_results_path = f'{results_path}/BASECALLS/{sample_name}'
    alignment_results_path = f'{results_path}/ALIGNMENTS/{sample_name}'
    variant_results_path = f'{results_path}/VARIANT_CALLS/{sample_name}'
    modification_results_path = f'{results_path}/MODIFICATION_COUNTS/{sample_name}'
    quality_control_results_path = f'{results_path}/QUALITY_CONTROL/{sample_name}'
    
    if not os.path.exists(basecall_results_path):
        os.makedirs(basecall_results_path)

    if not os.path.exists(alignment_results_path):
        os.makedirs(alignment_results_path)

    if not os.path.exists(variant_results_path):
        os.makedirs(variant_results_path)
        
    if not os.path.exists(modification_results_path):
        os.makedirs(modification_results_path)
        
    if not os.path.exists(quality_control_results_path):
        os.makedirs(quality_control_results_path)


    file_dict = {
        
        # BASECALL MODULE OUTFILES
        
        "basecall_file.bam": f'{basecall_results_path}/{sample_name}_basecall_with_{mod_name}_mods.bam',
        "summary_file.txt": f'{basecall_results_path}/{sample_name}_basecall_summary.txt',
        
        # ALIGNMENT MODULE OUTFILES
        
        "alignment_file.bam": f'{alignment_results_path}/{sample_name}_{ref_name}_alignment_with_{mod_name}_mods.bam',
        "filtered_alignment_file.bam": f'{alignment_results_path}/{sample_name}_{ref_name}_filtered_alignment_with_{mod_name}_mods.bam',
        
        # QUALITY CONTROL MODULE OUTFILES
        
        "alignment_pycoQC_file.html": f'{quality_control_results_path}/{sample_name}_{ref_name}_alignment_with_{mod_name}_mods_pycoQC.html',
        "filtered_alignment_pycoQC_file.html": f'{quality_control_results_path}/{sample_name}_{ref_name}_filtered_alignment_with_{mod_name}_mods_pycoQC.html',
        "alignment_qualimap_dir": f'{quality_control_results_path}/{sample_name}_{ref_name}_alignment_qualimap',
        "filtered_alignment_qualimap_dir": f'{quality_control_results_path}/{sample_name}_{ref_name}_filtered_alignment_qualimap',
        "merged_basecall_pycoQC": f"{results_path}/QUALITY_CONTROL/merged_basecall_pycoQC.html",
        "merged_alignment_pycoQC": f"{results_path}/QUALITY_CONTROL/merged_{ref_name}_alignment_pycoQC.html",
        
        # VARIANT MODULE OUTFILES
        
        "variant_call_dir": f'{variant_results_path}/{sample_name}_{ref_name}_filtered_alignment_variant_calls',
        "phased_alignment_file.bam": f'{variant_results_path}/{sample_name}_{ref_name}_phased_filtered_alignment_with_{mod_name}_mods.bam',
        
        # MODIFICATION COUNT MODULE OUTFILES
        
        "modkit_file.bed": f'{modification_results_path}/{sample_name}_{ref_name}_filtered_alignment_modkit.bed',
        "phased_modkit_file_dir": f'{modification_results_path}/{sample_name}_{ref_name}_phased_filtered_alignment_modkit'
        }


    ###############
    # BASECALLING #
    ###############


    gwf.target_from_template(
        name = f"basecalling_{sample_name}", 
        template = bctp.basecall(
            pod5_dir = f'{data_path}/POD5S/{sample_name}_pod5s/', 
            basecall_outfile = file_dict["basecall_file.bam"], 
            dorado_model = parameters["dorado_model"], 
            remora_model = parameters["remora_model"]))

    gwf.target_from_template(
        name = f"summary_{sample_name}", 
        template = bctp.create_summary_file(
            basecall_file = file_dict["basecall_file.bam"], 
            summary_file = file_dict["summary_file.txt"]))
    
    summary_file_list_pycoQC.append(file_dict["summary_file.txt"])


    ##############
    # ALIGNMENTS #
    ##############

    gwf.target_from_template(
        name = f"{ref_name}_alignment_{sample_name}", 
        template = altp.align_reads(
            basecall_infile = file_dict["basecall_file.bam"], 
            reference = parameters[f'reference'], 
            alignment_outfile = file_dict["alignment_file.bam"]))
    
    alignment_file_list_pycoQC.append(file_dict["alignment_file.bam"])
    
    gwf.target_from_template(
        name = f'filter_{ref_name}_alignment_{sample_name}', 
        template = altp.filter_alignment(
            alignment_infile = file_dict["alignment_file.bam"], 
            filtered_alignment_outfile = file_dict["filtered_alignment_file.bam"]))
    
    
    ###################
    # QUALITY CONTROL #
    ###################
    
    gwf.target_from_template(
        name = f"{ref_name}_alignment_qualimap_{sample_name}", 
        template = qctp.run_qualimap_on_alignment(
            alignment_infile = file_dict["alignment_file.bam"], 
            outfile_dir = file_dict["alignment_qualimap_dir"]))
    
    gwf.target_from_template(
        name = f"{ref_name}_filtered_alignment_qualimap_{sample_name}", 
        template = qctp.run_qualimap_on_alignment(
            alignment_infile = file_dict["filtered_alignment_file.bam"], 
            outfile_dir = file_dict["filtered_alignment_qualimap_dir"]))
    
    gwf.target_from_template(
        name = f"{ref_name}_pycoQC_{sample_name}", 
        template = qctp.run_pycoQC_on_alignment(
            summary_file = file_dict["summary_file.txt"], 
            alignment_file = file_dict["alignment_file.bam"], 
            QC_outfile = file_dict["alignment_pycoQC_file.html"], 
            QC_report_title = f'{sample_name}_alignment_with_{mod_name}_mods_pycoQC'))
    
    gwf.target_from_template(
        name = f"{ref_name}_filtered_alignment_pycoQC_{sample_name}", 
        template = qctp.run_pycoQC_on_alignment(
            summary_file = file_dict["summary_file.txt"], 
            alignment_file = file_dict["filtered_alignment_file.bam"], 
            QC_outfile = file_dict["filtered_alignment_pycoQC_file.html"], 
            QC_report_title = f'{sample_name}_filtered_alignment_with_{mod_name}_mods_pycoQC'))
    
    ###############################
    # VARIANT CALLING AND PHASING #
    ###############################
    
    gwf.target_from_template(
        name = f"call_SNPs_for_filtered_{ref_name}_alignment_{sample_name}", 
        template = vctp.call_SNPs(
            alignment_file = file_dict["filtered_alignment_file.bam"], 
            vcf_file_dir = file_dict["variant_call_dir"], 
            vcf_outfile = f'{file_dict["variant_call_dir"]}/merge_output.vcf.gz', 
            reference = parameters[f'reference'], 
            clair3_model_path = parameters[f'clair3_model']))
    
    gwf.target_from_template(
        name = f"phasing_SNPs_for_filtered_{ref_name}_alignment_{sample_name}", 
        template = vctp.phasing_whatshap(
            alignment_file = file_dict["filtered_alignment_file.bam"], 
            variant_file_dir = file_dict["variant_call_dir"], 
            reference = parameters[f'reference']))
    
    gwf.target_from_template(
        name = f"haplotagging_bam_for_filtered_{ref_name}_alignment_{sample_name}", 
        template = vctp.haplotagging_whatshap(
            alignment_file = file_dict["filtered_alignment_file.bam"], 
            reference = parameters[f'reference'], 
            variant_file_dir = file_dict["variant_call_dir"], 
            phased_bam_file = file_dict["phased_alignment_file.bam"]))
    
    ######################
    # MODIFICATION CALLS #
    ######################
    
    gwf.target_from_template(
        name = f"run_modkit_on_filtered_{ref_name}_alignment_{sample_name}", 
        template = mtp.run_modkit(
            alignment_file = file_dict["filtered_alignment_file.bam"], 
            modkit_outfile = file_dict["modkit_file.bed"], 
            reference = parameters[f'reference']))
    
    gwf.target_from_template(
        name = f"run_modkit_on_phased_filtered_{ref_name}_alignment_{sample_name}", 
        template = mtp.phasing_modkit(
            alignment_file = file_dict["phased_alignment_file.bam"], 
            reference = parameters[f'{ref_name}_reference'], 
            prefix = f'{sample_name}_phased_{ref_name}_modkit',
            outfile_dir = file_dict["phased_modkit_file_dir"]))
    
####################
# Merging QC files #
####################

gwf.target_from_template(
    name = f"merge_pycoQC_{ref_name}", 
    template = qctp.merge_pycoQCs(
        infile_list = summary_file_list_pycoQC, 
        QC_outfile = f"merged_basecall_{mod_name}_pycoQC.html", 
        QC_report_title = "merged_basecall_pycoQC"))

gwf.target_from_template(
    name = f"merge_alignment_pycoQC_{ref_name}", 
    template = qctp.merge_alignment_pycoQCs(
        summary_file_list = summary_file_list_pycoQC, 
        alignment_file_list = alignment_file_list_pycoQC, 
        QC_outfile = f"merged_{ref_name}_alignment_{mod_name}_pycoQC.html", 
        QC_report_title = "merged_alignment_pycoQC"))
