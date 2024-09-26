#def star_inputs(wildcards):
#    """Returns fastq inputs for star."""

#    base_path = "results/fastq/trimmed/{sample}.{lane}.{{pair}}.fastq.gz".format(
#        sample=wildcards.sample, lane=wildcards.lane
#    )
#    pairs = ["R1", "R2"] if is_paired else ["R1"]
#    print(expand(base_path, pair=pairs))
#    return expand(base_path, pair=pairs)



def star_extra(star_config):
    """Returns extra arguments for STAR based on config."""

    extra = star_config.get("extra", "")

    # Add readgroup information.
    extra_args = "--outSAMattrRGline " + star_config["readgroup"]

    # Add NM SAM attribute (required for PDX pipeline).
    if "--outSamAttributes" not in extra:
        extra_args += " --outSAMattributes NH HI AS nM NM"

    # Add any extra args passed by user.
    if extra:
        extra_args += " " + extra

    return extra_args


# create index in case missing
rule star_index:
    input:
        fasta="resources/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        directory(
            "resources/references/star_Homo_sapiens.GRCh38.dna.primary_assembly_overhang50"
        ),
    threads: 40
    params:
        extra="--sjdbGTFfile resources/references/Homo_sapiens.GRCh38.112.gtf --sjdbOverhang 50",
    log:
        "results/logs/star_index/hg38_human_overhang50.log",
    wrapper:
        "v4.5.0/bio/star/index"


rule star_align_pe:
    input:
        fq1="results/fastq/trimmed/{sample}.{lane}.R1.fastq.gz",
        fq2="results/fastq/trimmed/{sample}.{lane}.R2.fastq.gz",
        idx=config["star"]["index"],
    output: 
        aln="results/bam/star/{sample}.{lane}/Aligned.out.bam",            
        #log="results/logs/pe/{sample}.{lane}/Log.out",
        #sj="results/bam/star/pe/{sample}.{lane}/SJ.out.tab",
        #mapped=["results/bam/star/pe/{sample}.{lane}/unmapped.R1.fastq.gz","star/pe/{sample}.{lane}/unmapped.R2.fastq.gz"]
    log:
        "results/logs/star_align/{sample}.{lane}.log",
    params:
        extra="--outSAMtype BAM Unsorted",
    resources:
        memory=30,
    threads: config["star"]["threads"]
    wrapper:
        "v4.5.0/bio/star/align"


rule sambamba_sort:
    input:
        "results/bam/star/{sample}.{lane}/Aligned.out.bam"
    output:
        "results/bam/sorted/{sample}.{lane}.bam"
    params:
        config["sambamba_sort"]["extra"]
    threads:
        config["sambamba_sort"]["threads"]
    wrapper:
        "0.63.0/bio/sambamba/sort"


def merge_inputs(wildcards):
    lanes = get_sample_lanes(wildcards.sample)

    file_paths = [
        "results/bam/sorted/{}.{}.bam".format(wildcards.sample, lane) for lane in lanes
    ]

    return file_paths


rule samtools_merge:
    input:
        merge_inputs,
    output:
        "results/bam/final/{sample}.bam",
    params:
        config["samtools_merge"]["extra"],
    threads: config["samtools_merge"]["threads"]
    log:
        "results/logs/samtools_merge/{sample}.log",
    wrapper:
        "v4.5.0/bio/samtools/merge"


rule samtools_index:
    input:
        "results/bam/final/{sample}.bam",
    output:
        "results/bam/final/{sample}.bam.bai",
    log:
        "results/logs/samtools_index/{sample}.log",
    wrapper:
        "v4.5.0/bio/samtools/index"