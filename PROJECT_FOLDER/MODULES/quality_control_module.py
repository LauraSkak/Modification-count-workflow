from gwf import AnonymousTarget


def run_qualimap_on_alignment(alignment_infile, outfile_dir):
    """
    This creates a folder contain data and a quality control report for the alignment file using the qualimap software.
    """
    inputs = [alignment_infile]
    outputs = [f'{outfile_dir}/qualimapReport.html']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"36gb", "cores": 36}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate qualimap

    echo "Job ID: $SLURM_JOB_ID\n"
    
    echo "qualimap bamqc \
        -bam {infile} \
        -outdir {out_dir} \
        -nt 36 \
        --java-mem-size=32G\n"
    
    qualimap bamqc \
        -bam {infile} \
        -outdir {out_dir} \
        -nt 36 \
        --java-mem-size=32G
    
    '''.format(infile=alignment_infile, out_dir=outfile_dir)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_pycoQC_on_alignment(summary_file, alignment_file, QC_outfile, QC_report_title):
    """
    Takes an aligned .bam file produced by dorado aligner and a summary file produced by dorado summary to produce a quality control report using pycoQC software.
    """
    inputs = [summary_file, alignment_file]
    outputs = [QC_outfile]
    options = {"walltime":"05:00:00","account":"sexChromosomes", "memory":"50gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID\n"
    
    conda activate pycoQC
    
    echo "pycoQC --summary_file {sum_file} --bam_file {bam_file} --report_title {title} --html_outfile {outfile}\n"
    
    pycoQC \
        --summary_file {sum_file} \
        --bam_file {bam_file} \
        --report_title {title} \
        --html_outfile {outfile}
    
    '''.format(sum_file = summary_file, bam_file = alignment_file, outfile=QC_outfile, title = QC_report_title)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_pycoQCs(infile_list, QC_outfile, QC_report_title):
    """
    Takes a modbam file produced by the dorado basecaller and produces a quality control report using pycoQC software.
    """
    inputs = infile_list
    outputs = [QC_outfile]
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"200gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pycoQC

    echo "Job ID: $SLURM_JOB_ID\n"
    
    echo "pycoQC --summary_file {sum_files} --report_title {title} --html_outfile {outfile}\n"
    
    pycoQC \
        --summary_file {sum_files} \
        --report_title {title} \
        --html_outfile {outfile}
    
    '''.format(sum_files=" ".join(infile_list), outfile=QC_outfile, title = QC_report_title)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_alignment_pycoQCs(summary_file_list, alignment_file_list, QC_outfile, QC_report_title):
    """
    Takes a modbam file produced by the dorado basecaller and produces a quality control report using pycoQC software.
    """
    inputs = summary_file_list + alignment_file_list
    outputs = [QC_outfile]
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"500gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pycoQC

    echo "Job ID: $SLURM_JOB_ID\n"
    
    echo "pycoQC --summary_file {sum_files} --report_title {title} --html_outfile {outfile}\n"
    
    pycoQC \
        --summary_file {sum_files} \
        --bam_file {bam_files} \
        --report_title {title} \
        --html_outfile {outfile}
    
    '''.format(sum_files=" ".join(summary_file_list), bam_files=" ".join(alignment_file_list), outfile=QC_outfile, title = QC_report_title)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)