

def basecall(pod5_dir, basecall_outfile, dorado_model, remora_model):
    """
    Takes a directory containing pod5 files for a single sample and basecalls with modifications epigenetic modifications.
    """
    inputs = []
    outputs = [basecall_outfile]
    options = {"walltime":"48:00:00",
               "account":"sexChromosomes", 
               "memory":"25gb", 
               "gres":"gpu:1", 
               "queue":"gpu"}

    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "dorado basecaller {dorado_model} {pod5_dir} --modified-bases-models {remora_model} > {outfile}"

    dorado basecaller \
        {dorado_model} \
        {pod5_dir} \
        --modified-bases-models {remora_model} \
            > {outfile}
    
    '''.format(dorado_model=dorado_model, pod5_dir=pod5_dir, remora_model=remora_model, outfile=basecall_outfile)
    
    return inputs, outputs, options, spec


def run_fastQC_on_basecalls(basecall_infile, outfile_dir, outfile_prefix):
    """
    Takes a modbam file produced by the dorado basecaller and produces a quality control report using fastqc software.
    """
    inputs = [basecall_infile]
    outputs = [f'{outfile_dir}/{outfile_prefix}_fastqc.html']
    options = {"walltime":"05:00:00",
               "account":"sexChromosomes", 
               "memory":"10gb", 
               "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID\n"
    echo "samtools bam2fq --threads 16 {input} > {outfile_dir}/{outfile_prefix}.fastq\n"
    
    samtools bam2fq --threads 16 {input} > {outfile_dir}/{outfile_prefix}.fastq

    echo "fastqc --threads 16 -o {outfile_dir} {outfile_dir}/{outfile_prefix}.fastq\n"
    
    fastqc \
        --threads 16 \
        -o {outfile_dir} \
        {outfile_dir}/{outfile_prefix}.fastq

    rm {outfile_dir}/{outfile_prefix}.fastq
    
    '''.format(input=basecall_infile, outfile_dir=outfile_dir, outfile_prefix = outfile_prefix)
    
    return inputs, outputs, options, spec


def run_pycoQC_on_basecalls(basecall_infile, sum_file, QC_outfile, QC_report_title):
    """
    Takes a modbam file produced by the dorado basecaller and produces a quality control report using pycoQC software.
    """
    inputs = [basecall_infile]
    outputs = [QC_outfile]
    options = {"walltime":"01:00:00","account":"sexChromosomes", "memory":"50gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID\n"
    
    echo "dorado summary {infile} > {sum_file}\n"
    
    dorado summary {infile} > {sum_file}
    
    echo "pycoQC -f {sum_file} --report_title {title} -o {outfile}\n"
    
    pycoQC \
        -f {sum_file} \
        --report_title {title} \
        -o {outfile}
    
    '''.format(infile=basecall_infile, sum_file = sum_file, outfile=QC_outfile, title = QC_report_title)
    
    return inputs, outputs, options, spec