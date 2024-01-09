
def fai_index_reference(reference, fai_reference):
    """
    FIXME
    """
    inputs = [reference]
    outputs = [fai_reference]
    options = {"walltime":"01:00:00","account":"sexChromosomes", "memory":"16gb"}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "samtools faidx {ref} -o {fai_ref}"
    
    samtools faidx \
        {ref}  \
        -o {fai_ref} 
    
    '''.format(ref = reference, fai_ref = fai_reference)
    
    return inputs, outputs, options, spec


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
    conda activate epigen
    
    echo "Job ID: $SLURM_JOB_ID\n"
    echo "dorado aligner -t 32 {ref} {infile} | samtools sort > {outfile}\n"
    
    dorado aligner \
        --threads 32 \
        {ref} \
        {infile} \
        | samtools sort \
            > {outfile}
    
    echo "samtools index -@ 32 {outfile}\n"

    samtools index \
        -@ 32 \
        {outfile}
    
    '''.format(ref=reference, infile=basecall_infile, outfile=alignment_outfile)
    
    return inputs, outputs, options, spec


def filter_alignment(alignment_infile, filtered_alignment_outfile):
    """
    FIXME: Should the duplicated reads be filtered. It should them be --exclude-flags 1284
    """
    inputs = [alignment_infile]
    outputs = [filtered_alignment_outfile, f'{filtered_alignment_outfile}.bai']
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

    echo "Job ID: $SLURM_JOB_ID\n"

    echo "samtools view --threads 16 --bam --with-header --min-qlen 500 --exclude-flags 260 --expr "avg(qual) >= 10 && qlen >= 500" --min-MQ 10 {infile} > {outfile}"
    
    samtools view \
        --threads 16 \
        --bam --with-header \
        --min-qlen 500 \
        --exclude-flags 260 \
        --expr "avg(qual) >= 10 && qlen >= 500" \
        --min-MQ 10 \
        {infile} \
    > {outfile}
    
    echo "samtools index {outfile}\n"
    
    samtools index {outfile}
    
    '''.format(infile=alignment_infile, outfile=filtered_alignment_outfile)
    
    return inputs, outputs, options, spec

def run_qualimap_on_alignment(alignment_infile, outfile_dir):
    """
    FIXME
    """
    inputs = [alignment_infile]
    outputs = [f'{outfile_dir}/qualimapReport.html']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"36gb", "cores": 36}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate epigen

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
    
    return inputs, outputs, options, spec