# Methylation-count-for-ONT-long-read-sequencing-workflow
This workflow takes Oxford nanopore long-read sequencing data from raw signal to a csv with modification count data for each relevant position in a reference genome. Software packages such as Dorado, Clair3, whatshap and modkit is used.

All software used in the workflow will have an indication of software version. This doesn't not mean that newer versions of the software will not work the same. It just means that a workflow using these exact version has been proven to work. 

## The workflow

This sections explain each command used in the workflow.

![Workflow flowchart]("Workflow flowchart.png")

Necessary parameters include:

- **Dorado basecalling model**
- **Remora modification model** 
- **Reference genome**

**Dorado model** and **Remora model** is chosen based on how the data was produced. Look and [DNA models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) to see available models.


### The basecalling module

In this module we have to commands necessary for doing base- and modificationcalling and subsequent quality control.

#### Part 1: Basecalling

```shell
dorado basecaller \
    {dorado_model} \
    {pod5_dir} \
    --modified-bases-models {remora_model} \
        > {outfile}
```

This command uses the [Dorado software (version: 0.3.4+5f5cd02)](https://github.com/nanoporetech/dorado) to base- and modificationcall the raw .pod5 signal data to produce a .modbam file containing all basecalled reads with modification tags. 

- **dorado_model**: The full or relative path to the dorado DNA basecalling model.
- **pod5_dir**: The full or relative path to the folder containing all .pod5 associated with the sample.
- **remora_model**: The full or relative path to the remora modificationcalling model.
- **outfile**: The full or relative path to the .modbam file that should contain read data.

#### Part 2: Quality control

```shell
samtools bam2fq --threads 16 {infile} > {outfile_dir}/{outfile_prefix}.fastq

fastqc \
    --threads 16 \
    -o {outfile_dir} \
    {outfile_dir}/{outfile_prefix}.fastq

rm {outfile_dir}/{outfile_prefix}.fastq
```

This command uses the [Samtools software (version: 1.17)](https://github.com/samtools/samtools) to convert the .modbam file produced by dorado to a .fastq file. The .fastq file can then be used with the [fastQC software (v0.12.1)](https://github.com/s-andrews/FastQC) to create a quality control report. The fastQC software is primarily made for data produced for Illumina data, so some part of the report is irrelevant for Oxford nanopore data. The .fastq is removed after use. 16 cores are used here but no optimal running time analysis has been done, so it is not necessarily optimal.

- **infile**: The full or relative path to the .modbam file that contain read data produced in previous step.
- **outfile_dir**: The directory which should contain all files produced for the fastQC report.
- **outfile_prefix**: The desired fastQC report file name.

```shell
dorado summary {infile} > {sum_file}
    
pycoQC \
    -f {sum_file} \
    --report_title {title} \
    -o {outfile}

rm {sum_file}
```

Another option is the [pycoQC software (v2.5.0.3)](https://github.com/a-slide/pycoQC), which is a quality control software made for read data produced by ONTs old basecaller [Guppy](https://timkahlke.github.io/LongRead_tutorials/BS_G.html). Even though the software is not made for results produced by the dorado basecaller it seems to function properly.
First, a sequencing summary file is produced using the [Dorado software (version: 0.3.4+5f5cd02)](https://github.com/nanoporetech/dorado). This summary file is then used by the pycoQC software to produce a .html file with the quality control rapport.

- **infile**: The full or relative path to the .modbam file that contain read data produced in previous step.
- **sum_file**: The full or relative path to the intermediary summary file.
- **title**: The desired title of the quality control report.
- **outfile**: The full or relative path to the .html file containing the quality control report.


### The alignment module

After the raw data have been basecalled and a .modbam file with modification data has been produced, the sequencing reads can be aligned to a reference genome, whereafter the alignment can be filtered.  

#### Part 1: Reference genome preperation

```shell
samtools faidx \
        {ref}  \
        -o {fai_ref} 
```

A index file first has to be created for the reference genome. This can be done using the [Samtools software (version: 1.17)](https://github.com/samtools/samtools). 

- **ref**: The full or relative path to the .fasta file for the chosen reference genome.
- **fai_ref**: The full or relative path to the reference genome index file.

#### Part 2: Aligning the sequencing reads to the reference genome

```shell
dorado aligner \
    --threads 32 \
    {ref} \
    {infile} \
    | samtools sort \
        > {outfile}

samtools index \
    -@ 32 \
    {outfile}
```

The reads are first aligned to the reference using the dorado aligner which uses the ont-map method in minimap2 while also retaining modification data. 32 cores has been chosen as the thread count in this instance, but no running time analysis has been done to show this being the optimal choice. An index file is made to make viewing of resulting .bam file in [IGV](https://igv.org/) possible.

- **ref**: The full or relative path to the .fasta file for the chosen reference genome with an associated .fai fasta index file.
- **infile**: The full or relative path to the .modbam file containing read and modification data.
- **outfile**: The full or relative path to the desired location of the resulting .bam file.

#### Part 3: Filtering the alignment

```shell
samtools view \
    --threads 16 \
    --bam --with-header \
    --min-qlen 500 \
    --min-MQ 10 \
    --exclude-flags 260 \
    --expr "avg(qual) >= 10 && qlen >= 500" \
    {infile} \
> {outfile}

echo "samtools index {outfile}\n"

samtools index {outfile}
```

The alignment file produced by dorado aligner is filtered using [Samtools software (version: 1.17)](https://github.com/samtools/samtools). This filtering criteria used in this instance is; 

- A minimum query/mapped read length of 500 bp (--min-qlen 500) (remove????)
- A average mapping quality and minimum query/mapped read length of 500 bp (--expr "avg(qual) >= 10 && qlen >= 500")
- A minimum mapping quality of 10 (--min-MQ 10).
- only mapped and primary alignment is kept (and duplicates are exluded???)

NOTE TO SELF!: The minimum query length is set twice... is this a problem? remove one 
NOTE TO SELF!: The --exlude-flags tag value should probably be 1284, not 260 since 260 does not remove duplicates.

Filtering values are chosen based on the basecalling quality control and on how the data is later used for analysis.

- **infile**: The full or relative path to the .bam file produced by dorado aligner.
- **outfile**: The full or relative path to the desired location and name of the resulting filtered .bam file.

#### Part 4: Quality control

```shell
qualimap bamqc \
        -bam {infile} \
        -outdir {out_dir} \
        -nt 36 \
        --java-mem-size=32G
```

The [qualimap (v.2.2.2)](http://qualimap.conesalab.org/) is used for quality control of the alignment file.

- **infile**: The full or relative path to the filtered or unfiltered .bam file produced by dorado aligner.
- **outdir**: The full or relative path to the desired directory where the resulting folder containing all quality control results will be put.


### The variant calling and phasing module

If you want the methylation count data to be phased, you first need to phase the alignment file. This is done by using the Clair3 software to call SNPs whereafter the SNPs are phased and the .bam file is haplotagged using the whatshap software.

#### Part 1: Calling variants

```shell
run_clair3.sh \
    --bam_fn={infile} \
    --ref_fn={ref} \
    --output={outdir} \
    --threads=32 \
    --platform="ont" \
    --model_path={model_path} \
    --include_all_ctgs
```

The variants for all contigs/chromosomes are called using the [clair3 software (v1.0.4)](https://github.com/HKU-BAL/Clair3). A SNP calling model is chosen based on matching parameters to the dorado basecalling models. Available models can be found and downloaded with the [Rerio software](https://github.com/nanoporetech/rerio?tab=readme-ov-file#clair3-models).

- **infile**: The full or relative path to the filtered or unfiltered .bam file produced by dorado aligner.
- **ref**: The full or relative path to the reference genome fasta file.
- **outdir**: The full or relative path to the desired directory where the folder containing all intermediary and results files.
- **model_path**: The full or relative path to the clair3 model matching the same parameters chosen for the dorado basecalling model.

#### Part 2: Phasing reads

```shell
whatshap phase \
    --indels \
    --ignore-read-groups \
    --reference {ref} \
    --output {outfile} \
    {vcf_infile} {bam_infile}

tabix -p vcf {outfile}
```

All variants found, both SNPs and indels, are phased using the [whatshap software (version: 1.7)](https://github.com/whatshap/whatshap). All reads in .bam file come from the same sample, so any read-grouping is ignored. Lastly the phased vcf file is indexed using tabix.

- **vcf_infile**: The full or relative path to the merge_output.vcf.gz file found in the variant calling output directory.
- **bam_infile**: The full or relative path to the filtered or unfiltered .bam file produced by dorado aligner.
- **ref**: The full or relative path to the reference genome fasta file.
- **outfile**: The full or relative path to the desired output .vcf file.

#### Part 3: Haplotagging alignment file

```shell
whatshap haplotag \
    --out-threads 32 \
    --reference {ref} \
    --ignore-read-groups \
    --skip-missing-contigs \
    --tag-supplementary \
    {vcf_infile} {bam_infile} \
| samtools view \
    --threads 32 \
    --bam --with-header \
> {outfile}

samtools index -@ 32 {outfile}

whatshap stats --gtf {haplo_block_outfile} {vcf_infile}
```

The .bam file is haplotagged using [whatshap software (version: 1.7)](https://github.com/whatshap/whatshap) and the result is a new phased .bam file. It is made sure that all tags, such as modification tags are retained in the resulting phased .bam file. An index is created for the phased .bam file to make it possible to review in IGV. Lastly a .gtf file is created for the .bam file containing meta data on the haplotype blocks.

- **vcf_infile**: The full or relative path to the phased .vcf file produced by whatshap.
- **bam_infile**: The full or relative path to the filtered or unfiltered .bam file produced by dorado aligner.
- **ref**: The full or relative path to the reference genome fasta file.
- **outfile**: The full or relative path to the desired output .vcf file.
- **haplo_block_outfile**: The full or relative path to the desired output .gtf file.


### The modification call module

The last step in the module is precuring the methylation count data. This is done using the modkit software. 

#### For unphased data

```shell
modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --only-tabs \
        --threads 16
```

You use the [modkit software (version: 0.1.13)](https://github.com/nanoporetech/modkit) to create your modification count data. Modkit is the software recommended by ONT for methylation calling. If you have called more than one modification type, like both 5mC and 5hmC, when one line is created for each modification type and also for each strand (+ or -). To merge either strands or the modification types you can add the flags: 

    --combine-mods  (to merge the modification type)
    --combine-strands (to merge the data for each strand)

You can also use the --cpg flag if you only want data form reference CpG sites.

- **infile**: The full or relative path to the filtered .bam file produced by dorado aligner.
- **ref**: The full or relative path to the reference genome fasta file.
- **outfile**: The full or relative path to the desired output .csv file.

#### For phased data

```shell
modkit pileup \
        {infile} \
        {out_dir} \
        --ref {ref} \
        --partition-tag HP \
        --prefix {prefix} \
        --bedgraph \
        --only-tabs \
        --threads 16
```

If your read data is phased additional flags should be used to make seperate .csv file for each haplotype.

- **infile**: The full or relative path to the filtered .bam file produced by dorado aligner.
- **ref**: The full or relative path to the reference genome fasta file.
- **outdir**: The full or relative path to the directory which all resulting .csv file should be put.
- **prefix**: The desired prefix for the .cvs files.

Again the additional flags mentioned before can be added to the command;

    --combine-mods  (to merge the modification type)
    --combine-strands (to merge the data for each strand)
    --cpg




