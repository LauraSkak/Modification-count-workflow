from gwf import AnonymousTarget


def fai_index_reference(reference, fai_reference):
    """
    This function make a fasta index file for the input reference genome.
    """
    inputs = [reference]
    outputs = [fai_reference]
    options = {"walltime":"01:00:00","account":"sexChromosomes", "memory":"16gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "samtools faidx {ref} -o {fai_ref}"
    
    samtools faidx \
        {ref}  \
        -o {fai_ref} 
    
    '''.format(ref = reference, fai_ref = fai_reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def align_reads(basecall_infile, reference, alignment_outfile):
    """
    Takes the modbam file produced with dorado basecaller and aligns it to the specified reference using map-ont in minimap2. It produces a sorted bam file and a index file.
    """
    inputs = [basecall_infile]
    outputs = [alignment_outfile, f'{alignment_outfile}.bai']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"100gb", "cores": 32}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONTmod
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "dorado aligner -t 32 {ref} {infile} | samtools sort > {outfile}\n"
    
    dorado aligner \
        --threads 32 \
        {ref} \
        {infile} \
        | samtools sort \
            > {outfile}
    
    conda activate samtools
    
    echo "samtools index -@ 32 {outfile}\n"

    samtools index \
        -@ 32 \
        {outfile}
    
    '''.format(ref=reference, infile=basecall_infile, outfile=alignment_outfile)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    

def filter_alignment(alignment_infile, filtered_alignment_outfile):
    """
    This function filters the .bam file produced by dorado aligner. Its filters for reads that are a minimum for 500 bp in length, are the primary mapped read to the reference genome (no duplications), has a minimum mapping quality phred score of 10 and an average read quality phred score above 10.
    """
    inputs = [alignment_infile]
    outputs = [filtered_alignment_outfile, f'{filtered_alignment_outfile}.bai']
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools

    echo "Job ID: $SLURM_JOB_ID\n"
    
    samtools view \
        --threads 16 \
        --bam --with-header \
        --min-qlen 500 \
        --excl-flags 1284 \
        --expr "avg(qual) >= 10" \
        --min-MQ  10 \
        {infile} > {outfile}
    
    echo "samtools index {outfile}\n"
    
    samtools index {outfile}
    
    '''.format(infile=alignment_infile, outfile=filtered_alignment_outfile)
        
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#     echo "samtools view --threads 16 --bam --with-header --min-qlen 500 --excl-flags 1284 --expr "avg(qual) >= 10" --min-MQ  10 {infile} > {outfile}"