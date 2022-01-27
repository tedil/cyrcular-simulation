# TODO: either use --median_len and --sd_len params, or leave them to be determined from the trained error model (the latter is current behaviour)


rule merge_reads:
    input:
        fastqs=expand(
            "results/simulation/circular_reads/{ref}/{model}/{coverage}/{ref}.fastq",
            ref=list(regions["name"]),
            allow_missing=True,
        )
        + ["results/simulation/wgs_reads/whole_genome/{model}/whole_genome.fastq"],
    output:
        "results/simulation/merged_reads/{model}/{coverage}/reads.fastq.gz",
    shell:
        """cat {input.fastqs} | pigz -c > {output}"""


rule simulate_wgs_reads:
    input:
        reference=config["simulation"]["reference"]["path"],
        error_model="results/simulation/model/{model}/{model}.dir",
    output:
        out=directory("results/simulation/wgs_reads/whole_genome/{model}/"),
        error_profile="results/simulation/wgs_reads/whole_genome/{model}/whole_genome_aligned_error_profile",
        reads="results/simulation/wgs_reads/whole_genome/{model}/whole_genome.fastq",
    params:
        n_reads=lambda wc: int(
            np.ceil(
                config["simulation"]["wgs_noise"]["coverage"]
                * reference_size(wc)
                / median_readlength(wc)
            )
        ),
    conda:
        "../envs/nanosim.yaml"
    shell:
        """
        simulator.py genome --ref_g {input.reference} -c results/simulation/model/{wildcards.model} -o {output.out}/whole_genome -n {params.n_reads} -dna_type linear --fastq -t {threads} -b guppy
        cat {output.out}/whole_genome_aligned_reads.fastq {output.out}/whole_genome_unaligned_reads.fastq > {output.out}/whole_genome.fastq
        """


rule simulate_circular_reads_highest_coverage_only:
    input:
        reference="results/simulation/regions/{name}.fasta",
        error_model="results/simulation/model/{model}/{model}.dir",
    output:
        out=directory(
            "results/simulation/circular_reads/{name}/{model}/" + f"{highest_coverage}"
        ),
        error_profile="results/simulation/circular_reads/{name}/{model}/"
        + f"{highest_coverage}"
        + "/{name}_aligned_error_profile",
        reads="results/simulation/circular_reads/{name}/{model}/"
        + f"{highest_coverage}"
        + "/{name}.fastq",
    threads: 46
    conda:
        "../envs/nanosim.yaml"
    params:
        # TODO get actual median read length from error model directory
        n_reads=lambda wc: int(
            np.ceil(
                float(highest_coverage)
                * float(regions.loc[wc.name]["length"])
                / float(median_readlength(wc))
            )
        ),
    shell:
        """
        simulator.py genome --ref_g {input.reference} -c results/simulation/model/{wildcards.model} -o {output.out}/{wildcards.name} -n {params.n_reads} -dna_type circular  --fastq -t {threads} -b guppy
        cat {output.out}/{wildcards.name}_aligned_reads.fastq {output.out}/{wildcards.name}_unaligned_reads.fastq > {output.out}/{wildcards.name}.fastq
        """


rule simulate_circular_reads_by_subsampling:
    input:
        reference="results/simulation/regions/{name}.fasta",
        error_model="results/simulation/model/{model}/{model}.dir",
        reads="results/simulation/circular_reads/{name}/{model}/"
        + f"{highest_coverage}"
        + "/{name}.fastq",
    output:
        reads="results/simulation/circular_reads/{name}/{model}/{coverage}/{name}.fastq",
    threads: 1
    conda:
        "../envs/seqtk.yaml"
    params:
        # TODO get actual median read length from error model directory
        n_reads=lambda wc: int(
            np.ceil(
                float(wc.coverage)
                * float(regions.loc[wc.name]["length"])
                / float(median_readlength(wc))
            )
        ),
        seed=272977087,
    shell:
        """
        seqtk sample -2 -s {params.seed} {input.reads} {params.n_reads} > {output.reads}
        """


ruleorder: simulate_circular_reads_highest_coverage_only > simulate_circular_reads_by_subsampling


# TODO check for extension and unzip if needed
rule pipe_model_fastq:
    input:
        reads=config["simulation"]["error_model"]["sample"],
    output:
        pipe("tmp/{model}_sample.fastq"),
    shell:
        """
        cat {input.reads} > {output}
        """


rule build_error_model:
    input:
        reference=config["simulation"]["reference"]["path"],
        reads="tmp/{model}_sample.fastq",
    output:
        model=directory("results/simulation/model/{model}"),
        token="results/simulation/model/{model}/{model}.dir",
    threads: 48
    conda:
        "../envs/nanosim.yaml"
    shell:
        """
        read_analysis.py genome --read {input.reads} --ref_g {input.reference} -o {output.model} -t {threads} && touch {output.token}
        """


rule extract_region:
    output:
        "results/simulation/regions/{name}.fasta",
    params:
        ref=config["simulation"]["reference"]["path"],
        chrom=lambda wc: regions.loc[wc.name]["chrom"],
        start=lambda wc: regions.loc[wc.name]["start"],
        stop=lambda wc: regions.loc[wc.name]["stop"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/extract_region/{name}.log",
    shell:
        """samtools faidx {params.ref} {params.chrom}:{params.start}-{params.stop} > {output} 2> {log}"""