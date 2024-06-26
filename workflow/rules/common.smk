from snakemake.utils import validate
import numpy as np
import pandas as pd
import gzip
from functools import lru_cache

validate(config, schema="../schemas/config.schema.yaml")

regions = pd.read_csv(
    config["simulation"]["regions"],
    sep="\t",
    comment="#",
    dtype={"name": str, "chrom": str, "start": int, "stop": int, "length": int},
)

regions["length"] = regions["stop"] - regions["start"]
regions = regions.set_index("name", drop=False)
validate(regions, schema="../schemas/regions.schema.yaml")

highest_coverage = np.max(config["simulation"]["target_coverages"])


def simulation_target_params():
    model_name = config["simulation"]["error_model"]["nanopore"]["name"]
    models = [model_name]
    model_name_illumina = (
        config["simulation"]["error_model"].get("illumina", {}).get("name", None)
    )
    if model_name_illumina:
        models += [model_name_illumina]
    coverages = config["simulation"]["target_coverages"]
    return {
        model: {
            coverage: [
                f"results/simulation/merged_reads/{model}/{coverage}/reads.fastq.gz"
            ]
            if model != model_name_illumina
            else [
                f"results/simulation/merged_reads/{model}/{coverage}/reads_1.fastq.gz",
                f"results/simulation/merged_reads/{model}/{coverage}/reads_2.fastq.gz",
            ]
            for coverage in coverages
        }
        for model in models
    }


def simulation_targets():
    params = simulation_target_params()
    return [
        f
        for model, info in params.items()
        for _coverage, fq in info.items()
        for f in fq
    ]


def make_samplesheet():
    params = simulation_target_params()
    s = []
    for model, info in params.items():
        for coverage, fqs in info.items():
            group = sample = f"{model}_{coverage}"
            s.append([group, sample, "nanopore"])
    samples = pd.DataFrame(columns=["group", "sample", "platform"], data=s)
    print(samples)
    return samples


samples = make_samplesheet()


def make_unitsheet():
    params = simulation_target_params()
    print(params)
    u = []
    for model, info in params.items():
        for coverage, fqs in info.items():
            group = sample = f"{model}_{coverage}"
            if len(fqs) == 1:
                u.append([sample, "0", fqs[0], None])
            elif len(fqs) == 2:
                u.append([sample, "0", fqs[0], fqs[1]])
            else:
                raise ValueError(
                    f"Expected either single or paired-end data but got {len(fqs)} fastq files"
                )
    units = pd.DataFrame(columns=["sample", "unit", "fq1", "fq2"], data=u)
    print(units)
    return units


units = make_unitsheet()

global median_read_length
global reference_length

median_read_length = None
reference_length = None


def median_readlength(wildcards):
    global median_read_length
    if median_read_length is not None:
        return median_read_length
    # FIXME this should extract the median read length from the nanosim model
    sample_path = config["simulation"]["error_model"]["sample"]
    opener = gzip.open if sample_path.endswith(".gz") else open
    lengths = []
    with opener(sample_path, "r") as f:
        for i, line in enumerate(f.readlines()):
            if i % 4 == 1:
                l = len(line) - 1  # do not include newline char in length
                lengths.append(l)
    median_read_length = np.median(lengths)
    return median_read_length


@lru_cache(maxsize=2)
def reference_size(wildcards):
    global reference_length
    if reference_length is not None:
        return reference_length
    ref_path = config["simulation"]["reference"]["path"]
    opener = gzip.open if ref_path.endswith(".gz") else open
    description_contains_length = True
    length = 0
    with opener(ref_path, "r") as f:
        for line in f.readlines():
            if description_contains_length:
                line = line.strip()
                if line.startswith(">"):
                    if "LN:" in line:
                        part = next(
                            filter(lambda p: p.startswith("LN:"), line.split(" "))
                        )
                        if part:
                            l = int(part.lstrip("LN:"))
                            length += l
                        else:
                            raise ValueError("failed parsing length")
                    else:
                        description_contains_length = False
                else:
                    continue
            else:
                if line.startswith(">"):
                    continue
                length += len(line) - 1
    print(length)
    reference_length = length
    return reference_length


def all_input(wildcards):
    model_name = config["simulation"]["error_model"]["name"]
    inputs = ["results/calling/config/samples.tsv", "results/calling/config/units.tsv"]
    inputs += [f"results/simulation/model/{model_name}/{model_name}.dir"]
    inputs += simulation_targets()
    inputs += [f"results/validation/{model_name}.stats.tsv"]
    return inputs
