
from gwf import AnonymousTarget

def run_modkit(alignment_file, modkit_outfile, reference):
    """
    Takes an unphased .bam file and creates a methylation count table.
    """
    inputs = [alignment_file]
    outputs = [modkit_outfile]
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16"

    modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --only-tabs \
        --threads 16
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def phasing_modkit(alignment_file, reference, prefix, outfile_dir):
    """
    Takes a phased .bam file and creates a methylation count table for each haplotype.
    """
    inputs = [alignment_file, f'{alignment_file}.bai']
    outputs = [f'{outfile_dir}/{prefix}.bed']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

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
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_EPIC_modkit(alignment_file, modkit_outfile, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strand and modification types are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [modkit_outfile]
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16"

    modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --threads 16 \
        --combine-mods \
        --combine-strands \
        --only-tabs \
        --cpg
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
