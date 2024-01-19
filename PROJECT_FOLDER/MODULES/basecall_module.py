from gwf import AnonymousTarget

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
    conda activate ONTmod
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "dorado basecaller {dorado_model} {pod5_dir} --modified-bases-models {remora_model} > {outfile}"

    dorado basecaller \
        {dorado_model} \
        {pod5_dir} \
        --modified-bases-models {remora_model} \
            > {outfile}
    
    '''.format(dorado_model=dorado_model, pod5_dir=pod5_dir, remora_model=remora_model, outfile=basecall_outfile)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def create_summary_file(basecall_file, summary_file):
    """
    Creates a summary file of the nanopore read data. Is used for producing pycoQC reports.
    """
    inputs = [basecall_file]
    outputs = [summary_file]
    options = {"walltime":"05:00:00","account":"sexChromosomes", "memory":"50gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID\n"
    
    echo "dorado summary {infile} > {sum_file}\n"
    
    dorado summary {infile} > {sum_file}
    
    '''.format(infile = basecall_file, sum_file = summary_file)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)