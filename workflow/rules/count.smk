from os import path
import numpy as np

def feature_counts_extra(wildcards):
    extra = config["feature_counts"]["extra"]
    if is_paired:
        extra += " -p"
    return extra


rule feature_counts:
    input:
        bam="results/bam/final/{sample}.bam",
        bai="results/bam/final/{sample}.bam.bai",
        annotation=config["feature_counts"]["annotation"]
    output:
        counts="results/counts/per_sample/{sample}.txt",
        summary="results/qc/feature_counts/{sample}.txt"
    params:
        extra=feature_counts_extra
    threads:
        config["feature_counts"]["threads"]
    log:
        "results/logs/feature_counts/{sample}.txt"
    wrapper:
            "v4.5.0/bio/subread/featurecounts"


rule merge_counts:
    input:
        expand("results/counts/per_sample/{sample}.txt", sample=get_samples())
    output:
        "results/counts/merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,
                        index_col=list(range(6)))
            for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)
