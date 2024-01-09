
def run_modkit(alignment_file, modkit_outfile, reference):
    """
    FIXME
    Maybe use --bedgraph
    """
    inputs = [alignment_file]
    outputs = [modkit_outfile]
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16"

    modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --only-tabs \
        --threads 16
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return inputs, outputs, options, spec


def phasing_modkit(alignment_file, reference, prefix, outfile_dir):
    """
    FIXME 
    is used after phasing
    
    """
    inputs = [alignment_file, f'{alignment_file}.bai']
    outputs = [f'{outfile_dir}/{prefix}.bed']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {prefix}.bed --ref {ref} --partition-tag HP --prefix {prefix} --threads 16"

    modkit pileup \
        {infile} \
        {out_dir} \
        --ref {ref} \
        --partition-tag HP \
        --prefix {prefix} \
        --bedgraph \
        --only-tabs \
        --threads 16
    
    '''.format(ref = reference, infile = alignment_file, prefix=prefix, out_dir = outfile_dir)
    
    return inputs, outputs, options, spec

def run_EPIC_modkit(alignment_file, modkit_outfile, reference):
    """
    FIXME
    """
    inputs = [alignment_file]
    outputs = [modkit_outfile]
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16"

    modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --threads 16 \
        --combine-mods \
        --combine-strands \
        --cpg
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return inputs, outputs, options, spec