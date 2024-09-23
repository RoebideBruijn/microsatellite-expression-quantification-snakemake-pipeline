from os import path

def multiqc_inputs(wildcards):
    """Returns inputs for multiqc, which vary depending on whether pipeline
       is processing normal or PDX data and whether the data is paired."""

    inputs = [
        expand("results/qc/fastqc/{sample_lane}.{pair}_fastqc.html",
               sample_lane=get_samples_with_lane(),
               pair=["R1", "R2"] if is_paired else ["R1"]),
        expand("results/qc/cutadapt/{sample_lane}.txt",
               sample_lane=get_samples_with_lane()),
        expand("results/qc/samtools_stats/{sample}.txt", sample=get_samples())
    ]

    return [input_ for sub_inputs in inputs for input_ in sub_inputs]

rule multiqc:
    input:
        multiqc_inputs
    output:
        "results/qc/multiqc_report.html"
    params:
        config["multiqc"]["extra"]
    log:
        "resuls/logs/multiqc.log"
    wrapper:
        "v4.5.0/bio/multiqc"

rule fastqc:
    input:
        "results/fastq/trimmed/{sample}.{lane}.{pair}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.{lane}.{pair}_fastqc.html",
        zip="results/qc/fastqc/{sample}.{lane}.{pair}_fastqc.zip"
    params:
        config["fastqc"]["extra"]
    wrapper:
        "v4.5.0/bio/fastqc"

rule samtools_stats:
    input:
        "results/bam/final/{sample}.bam"
    output:
        "results/qc/samtools_stats/{sample}.txt"
    wrapper:
        "v4.5.0/bio/samtools/stats"