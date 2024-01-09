
def call_SNPs(alignment_file, vcf_file_dir, vcf_outfile, reference, clair3_model_path):
    """
    Takes an alignment file and creates a vcf file containing all SNPs found in the .bam file.
    
    TODO: Maybe i should add the --sample_name tag, so i can later merge the vcf files for all samples.
    """
    inputs = [alignment_file, f'{alignment_file}.bai', reference]
    outputs = [f'{vcf_outfile}.tbi', vcf_outfile]
    options = {"walltime":"24:00:00","account":"sexChromosomes", "memory":"100gb", "cores": 32}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate clair3

    echo "Job ID: $SLURM_JOB_ID\n"

    echo "run_clair3.sh --bam_fn={infile} --ref_fn={ref} --output={out_dir} --threads=32 --platform="ont" --model_path={model_path} --include_all_ctgs"
    
    run_clair3.sh \
        --bam_fn={infile} \
        --ref_fn={ref} \
        --output={out_dir} \
        --threads=32 \
        --platform="ont" \
        --model_path={model_path} \
        --include_all_ctgs
    
    '''.format(infile = alignment_file, ref=reference, model_path = clair3_model_path, out_dir = vcf_file_dir)
    
    return inputs, outputs, options, spec


def phasing_whatshap(alignment_file, variant_file_dir, reference):
    """
    Phases a vcf file based on a .bam file produced by dorado aligner and a .vcf file produced with clair3. The .vcf file is then indexed with tabix, so the resulting files should be a phased .vcf file and a .tbi index file.
    """
    inputs = [alignment_file, f'{variant_file_dir}/merge_output.vcf.gz']
    outputs = [f'{variant_file_dir}/phased.vcf.gz', f'{variant_file_dir}/phased.vcf.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 1}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate clair3

    echo "Job ID: $SLURM_JOB_ID"
    
    echo "whatshap phase --indels --ignore-read-groups --reference {ref} --output {phased_variant_outfile} {vcf_infile} {bam_infile}"

    whatshap phase \
        --indels \
        --ignore-read-groups \
        --reference {ref} \
        --output {phased_variant_outfile} \
        {vcf_infile} {bam_infile}

    echo "tabix -p vcf {phased_variant_outfile}"

    tabix -p vcf {phased_variant_outfile}
    
    '''.format(ref = reference, phased_variant_outfile = f'{variant_file_dir}/phased.vcf.gz', vcf_infile = f'{variant_file_dir}/merge_output.vcf.gz', bam_infile = alignment_file)
    
    return inputs, outputs, options, spec

def haplotagging_whatshap(alignment_file, reference, variant_file_dir, phased_bam_file):
    """
    Makes a new .bam file containing a haplotag for each read. A .gtf file containing data for continously phased segments are also created.
    """
    inputs = [alignment_file, f'{variant_file_dir}/phased.vcf.gz', f'{variant_file_dir}/phased.vcf.gz.tbi']
    outputs = [phased_bam_file, f'{variant_file_dir}/haplo_blocks.gft']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 32}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate clair3

    echo "Job ID: $SLURM_JOB_ID"
    
    echo "whatshap haplotag --out-threads 32 --reference {ref} --ignore-read-groups --skip-missing-contigs --tag-supplementary {phased_variant_outfile} {bam_infile} | samtools view --threads 32 --bam --with-header > {phased_bam_outfile}"

    whatshap haplotag \
        --out-threads 32 \
        --reference {ref} \
        --ignore-read-groups \
        --skip-missing-contigs \
        --tag-supplementary \
        {phased_variant_outfile} {bam_infile} \
    | samtools view \
        --threads 32 \
        --bam --with-header \
    > {phased_bam_outfile}

    echo "samtools index -@ 32 {phased_bam_outfile}"
    
    samtools index -@ 32 {phased_bam_outfile}

    echo "whatshap stats --gtf {haplo_block_outfile} {phased_variant_outfile}"
    
    whatshap stats --gtf {haplo_block_outfile} {phased_variant_outfile}
    
    '''.format(ref = reference, phased_variant_outfile = f'{variant_file_dir}/phased.vcf.gz', bam_infile = alignment_file, phased_bam_outfile = phased_bam_file, haplo_block_outfile = f'{variant_file_dir}/haplo_blocks.gft')
    
    return inputs, outputs, options, spec