from os import path


def input_path(wildcards):
    """Extracts input path from sample overview."""

    # Lookup sample for given lane/sample ids.
    subset = samples.query(
        "sample == {!r} and lane == {!r}".format(wildcards.sample, wildcards.lane)
    )
    # Extract file_path.
    fastq_col = "fastq1" if wildcards.pair == "R1" else "fastq2"
    file_path = subset.iloc[0][fastq_col]

    # Prepend local directory if given.
    input_dir = config["input"].get("dir", None)

    if input_dir is not None:
        file_path = path.join(input_dir, file_path)

    return file_path


rule copy_input:
    input:
        input_path,
    output:
        temp("results/fastq/raw/{sample}.{lane}.{pair}.fastq.gz"),
    resources:
        io=1,
    shell:
        "cp {input} {output}"