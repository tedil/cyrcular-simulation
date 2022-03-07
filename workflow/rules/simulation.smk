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
    conda:
        "../envs/pigz.yaml"
    log:
        "logs/merge_reads/{model}_{coverage}.log",
    shell:
        """cat {input.fastqs} | pigz -c > {output} 2> {log}"""


rule merge_reads_illumina:
    input:
        fq1=expand(
            "results/simulation/circular_reads/{ref}/illumina/{coverage}/{ref}_1.fastq.gz",
            ref=list(regions["name"]),
            allow_missing=True,
        )
        + [
            f'results/simulation/wgs_reads/whole_genome/illumina/{config["simulation"]["reference"]["name"]}_1.fastq.gz'
        ],
        fq2=expand(
            "results/simulation/circular_reads/{ref}/illumina/{coverage}/{ref}_2.fastq.gz",
            ref=list(regions["name"]),
            allow_missing=True,
        )
        + [
            f'results/simulation/wgs_reads/whole_genome/illumina/{config["simulation"]["reference"]["name"]}_2.fastq.gz',
        ],
    output:
        fq1="results/simulation/merged_reads/{model}/{coverage}/reads_1.fastq.gz",
        fq2="results/simulation/merged_reads/{model}/{coverage}/reads_2.fastq.gz",
    conda:
        "../envs/pigz.yaml"
    log:
        "logs/merge_reads/{model}_{coverage}.log",
    shell:
        """pigz -dc {input.fq1} | pigz -c > {output.fq1} 2> {log}
        pigz -dc {input.fq2} | pigz -c > {output.fq2} 2>> {log}
        """


rule simulate_circular_reads_illumina_highest_coverage_only:
    input:
        reference="results/simulation/regions/{name}.fasta",
    output:
        out=directory(
            f"results/simulation/circular_reads/{{name}}/illumina/{highest_coverage}"
        ),
        fq1=f"results/simulation/circular_reads/{{name}}/illumina/{highest_coverage}/{{name}}_1.fastq.gz",
        fq2=f"results/simulation/circular_reads/{{name}}/illumina/{highest_coverage}/{{name}}_2.fastq.gz",
    threads: 1
    conda:
        "../envs/readSimulator.yaml"
    log:
        "logs/simulate_circular_reads_illumina_highest_coverage_only/{name}_illumina.log",
    params:
        depth=highest_coverage,
        seed=272977087,
    shell:
        """
        readSimulator.py \
        --input {input.reference} \
        --simulator art \
        --simulator_path $(which art_illumina) \
        --outdir {output.out} \
        --iterations 10 \
        --readlen 100 \
        --depth {params.depth} \
        --opts '-m 300 -s 25 -rs {params.seed} -ss HS20 -na'
        """


rule simulate_circular_reads_illumina_by_subsampling:
    input:
        reference="results/simulation/regions/{name}.fasta",
        fq1=f"results/simulation/circular_reads/{{name}}/illumina/{highest_coverage}/{{name}}_1.fastq.gz",
        fq2=f"results/simulation/circular_reads/{{name}}/illumina/{highest_coverage}/{{name}}_2.fastq.gz",
    output:
        fq1="results/simulation/circular_reads/{name}/illumina/{coverage}/{name}_1.fastq.gz",
        fq2="results/simulation/circular_reads/{name}/illumina/{coverage}/{name}_2.fastq.gz",
    threads: 1
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/simulate_circular_reads_illumina_by_subsampling/{name}_illumina_{coverage}.log",
    params:
        fraction=lambda wc: float(wc.coverage) / float(highest_coverage),
        seed=272977087,
    shell:
        """
        seqtk sample -2 -s {params.seed} {input.fq1} {params.fraction} > {output.fq1} 2> {log}
        seqtk sample -2 -s {params.seed} {input.fq2} {params.fraction} > {output.fq2} 2>> {log}
        """


rule simulate_wgs_reads_illumina:
    input:
        reference=config["simulation"]["reference"]["path"],
    output:
        out=directory("results/simulation/wgs_reads/whole_genome/illumina/"),
        fq1=f'results/simulation/wgs_reads/whole_genome/illumina/{config["simulation"]["reference"]["name"]}_1.fastq.gz',
        fq2=f'results/simulation/wgs_reads/whole_genome/illumina/{config["simulation"]["reference"]["name"]}_2.fastq.gz',
    params:
        depth=float(config["simulation"]["wgs_noise"]["coverage"]),
        seed=272977087,
    log:
        "logs/simulate_wgs_reads/illumina.log",
    conda:
        "../envs/readSimulator.yaml"
    shell:
        """
        readSimulator.py \
        --input {input.reference} \
        --simulator art \
        --simulator_path $(which art_illumina) \
        --outdir {output.out} \
        --iterations 1 \
        --readlen 100 \
        --depth {params.depth} \
        --opts '-m 300 -s 25 -rs {params.seed} -ss HS20 -na'
        """


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
    log:
        "logs/simulate_wgs_reads/{model}.log",
    conda:
        "../envs/nanosim.yaml"
    shell:
        """
        simulator.py genome --ref_g {input.reference} -c results/simulation/model/{wildcards.model} -o {output.out}/whole_genome -n {params.n_reads} -dna_type linear --fastq -t {threads} -b guppy 2> {log}
        cat {output.out}/whole_genome_aligned_reads.fastq {output.out}/whole_genome_unaligned_reads.fastq > {output.out}/whole_genome.fastq 2>> {log}
        """


rule simulate_circular_reads_highest_coverage_only:
    input:
        reference="results/simulation/regions/{name}.fasta",
        error_model="results/simulation/model/{model}/{model}.dir",
    output:
        out=directory(
            f"results/simulation/circular_reads/{{name}}/{{model}}/{highest_coverage}"
        ),
        error_profile=f"results/simulation/circular_reads/{{name}}/{{model}}/{highest_coverage}/{{name}}_aligned_error_profile",
        reads=f"results/simulation/circular_reads/{{name}}/{{model}}/{highest_coverage}/{{name}}.fastq",
    threads: 46
    conda:
        "../envs/nanosim.yaml"
    log:
        "logs/simulate_circular_reads_highest_coverage_only/{name}_{model}.log",
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
        simulator.py genome --ref_g {input.reference} -c results/simulation/model/{wildcards.model} -o {output.out}/{wildcards.name} -n {params.n_reads} -dna_type circular  --fastq -t {threads} -b guppy 2> {log}
        cat {output.out}/{wildcards.name}_aligned_reads.fastq {output.out}/{wildcards.name}_unaligned_reads.fastq > {output.out}/{wildcards.name}.fastq 2>> {log}
        """


rule simulate_circular_reads_by_subsampling:
    input:
        reference="results/simulation/regions/{name}.fasta",
        error_model="results/simulation/model/{model}/{model}.dir",
        reads=f"results/simulation/circular_reads/{{name}}/{{model}}/{highest_coverage}/{{name}}.fastq",
    output:
        reads="results/simulation/circular_reads/{name}/{model}/{coverage}/{name}.fastq",
    threads: 1
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/simulate_circular_reads_by_subsampling/{name}_{model}_{coverage}.log",
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
        seqtk sample -2 -s {params.seed} {input.reads} {params.n_reads} > {output.reads} 2> {log}
        """


ruleorder: simulate_circular_reads_illumina_highest_coverage_only > simulate_circular_reads_highest_coverage_only > simulate_circular_reads_by_subsampling


# TODO check for extension and unzip if needed
rule pipe_nanopore_model_fastq:
    input:
        reads=config["simulation"]["error_model"]["nanopore"]["sample"],
    output:
        pipe("tmp/{model}_sample.fastq"),
    log:
        "logs/pipe_model_fastq/{model}.log",
    conda:
        "../envs/pigz.yaml"
    shell:
        """
        cat {input.reads} > {output} 2> {log}
        """


rule build_nanopore_error_model:
    input:
        reference=config["simulation"]["reference"]["path"],
        reads="tmp/{model}_sample.fastq",
    output:
        model=directory("results/simulation/model/{model}"),
        token="results/simulation/model/{model}/{model}.dir",
    threads: 48
    conda:
        "../envs/nanosim.yaml"
    log:
        "logs/build_error_model/{model}.log",
    shell:
        """
        read_analysis.py genome --read {input.reads} --ref_g {input.reference} -o {output.model} -t {threads} 2> {log} && touch {output.token}
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
