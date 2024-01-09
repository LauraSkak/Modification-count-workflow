# Methylation-count-for-ONT-long-read-sequencing-workflow
This workflow takes Oxford nanopore signal data from raw signal to a csv with modification data for each relevant position in a reference genome. Software like Dorado, Clair3, whatshap and modkit is used.

## The workflow

This sections explain each command used in the workflow.

![Workflow flowchart](Workflow flowchart.png)

Necessary parameters include:

- **Dorado basecalling model**
- **Remora modification model** 
- **Reference genome**

**Dorado model** and **Remora model** is chosen based on how the data was produced. Look and [DNA models](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) to see available models.

### The basecalling module

In this module we have to commands necessary for doing base- and modificationcalling and subsequent quality control.

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


```shell
samtools bam2fq --threads 16 {infile} > {outfile_dir}/{outfile_prefix}.fastq

fastqc \
    --threads 16 \
    -o {outfile_dir} \
    {outfile_dir}/{outfile_prefix}.fastq

rm {outfile_dir}/{outfile_prefix}.fastq
```

This command uses the [Samtools software (version: 1.17)](https://github.com/samtools/samtools) to convert the .modbam file produced by dorado to a .fastq file. The .fastq file can then be used with the [fastQC software (version: v0.12.1)](https://github.com/s-andrews/FastQC) to create a quality control report. The fastQC software is primarily made for data produced for Illumina data, so some part of the report is irrelevant for Oxford nanopore data. The .fastq is removed after use. 16 cores are used here but no optimal running time analysis has been done, so it is not necessarily optimal.

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

Another option is the [pycoQC software (version: v2.5.0.3)](https://github.com/a-slide/pycoQC), which is a quality control software made for read data produced by ONTs old basecaller [Guppy](https://timkahlke.github.io/LongRead_tutorials/BS_G.html). Even though the software is not made for results produced by the dorado basecaller it seems to function properly.
First, a sequencing summary file is produced using the [Dorado software (version: 0.3.4+5f5cd02)](https://github.com/nanoporetech/dorado). This summary file is then used by the pycoQC software to produce a .html file with the quality control rapport.

- **infile**: The full or relative path to the .modbam file that contain read data produced in previous step.
- **sum_file**: The full or relative path to the intermediary summary file.
- **title**: The desired title of the quality control report.
- **outfile**: The full or relative path to the .html file containing the quality control report.


### The alignment module

