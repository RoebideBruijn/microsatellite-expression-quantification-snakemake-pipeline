if is_paired:

    rule cutadapt:
        input:
            [
                "results/fastq/raw/{sample}.{lane}.R1.fastq.gz",
             "results/fastq/raw/{sample}.{lane}.R2.fastq.gz",
            ],
        output:
            fastq1=temp("results/fastq/trimmed/{sample}.{lane}.R1.fastq.gz"),
            fastq2=temp("results/fastq/trimmed/{sample}.{lane}.R2.fastq.gz"),
            qc="results/qc/cutadapt/{sample}.{lane}.txt",
        params:
            extra="-q 20 --minimum-length 20",
            adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        log:
            "results/logs/cutadapt/{sample}.{lane}.log",
        wrapper:
            "v4.5.0/bio/cutadapt/pe"

else:

    rule cutadapt:
        input:
            "results/fastq/raw/{sample}.{lane}.R1.fastq.gz",
        output:
            fastq=temp("results/fastq/trimmed/{sample}.{lane}.R1.fastq.gz"),
            qc="results/qc/cutadapt/{sample}.{lane}.txt",
        params:
            config["cutadapt_se"]["extra"],
        log:
            "results/logs/cutadapt/{sample}.{lane}.log",
        wrapper:
            "v4.5.0/bio/cutadapt/se"