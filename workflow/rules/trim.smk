rule cutadapt:
    input:
        [
        "results/fastq/raw/{sample}.{lane}.R1.fastq.gz",
        "results/fastq/raw/{sample}.{lane}.R2.fastq.gz",
        ],
    output:
        fastq1="results/fastq/trimmed/{sample}.{lane}.R1.fastq.gz",
        fastq2="results/fastq/trimmed/{sample}.{lane}.R2.fastq.gz",
        qc="results/qc/cutadapt/{sample}.{lane}.txt",
    params:
        extra="-q 20 --minimum-length 20",
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    log:
        "results/logs/cutadapt/{sample}.{lane}.log",
    wrapper:
        "v4.5.0/bio/cutadapt/pe"