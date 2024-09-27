# ======================================================================================================================
# Project: Project_ABC
# Script : build_baseline.smk
# Author : Peng Jia
# Date   :  2024/5/27
# Email  : pengjia@xjtu.edu.cn
# Description: Pipeline of baseline building of MSIsensor-pro
# ======================================================================================================================
configfile: "config.yaml"
import pandas as pd

print(config)
dir_work = config["dir_work"]
dir_work = dir_work if dir_work.endswith("/") else dir_work + "/"
baseline_bams = config["samples4baseline"]
baseline_samples = pd.read_table(f"{baseline_bams}",index_col=0)
samples = pd.read_table(f"{config['samples4test']}",index_col=0)
msisensor_pro = config["msisensor-pro"]
genome_version = config["genome_version"]
print(samples)
print(baseline_samples)
genome_reference = {config["genome_version"]: config["reference"]}
# print(dir_work + f"baselines/{genome_version}.baseline.file")

rule all:
    input:
        dir_work+ f"tumor_only_output.{genome_version}.merge.tsv"

        # dir_work + f"baselines/{genome_version}.baseline.file"


rule scan:
    input:
        lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        dir_work + "reference/{genome_version}.msisensor.scan.list"
    log:
        dir_work + "reference/{genome_version}.msisensor.scan.log"
    run:
        shell("{msisensor_pro} scan -d {input} -o {output} 2>{log} 1>{log}")

rule pro_pre:
    input:
        bam=lambda wildcards: baseline_samples.loc[wildcards.sample, "bam_path"],
        ms=dir_work + "reference/{genome_version}.msisensor.scan.list",
        ref=lambda wildcards: genome_reference[wildcards.genome_version],
    output:
        dir_work + "baselines/details/{sample}.{genome_version}.baseline.out"
    threads: 2
    run:
        shell("{msisensor_pro} pro -d {input.ms} -t {input.bam} -g {input.ref} -o {output}")

rule merge:
    input:
        expand(dir_work + "baselines/details/{sample}.{{genome_version}}.baseline.out",sample=baseline_samples.index)
    output:
        dir_work + "baselines/{genome_version}.baseline.samples.list"
    run:
        with open(f"{output}","w") as f:
            for i in input:
                f.write(f"{i.split('/')[-1]} \t{i}_all\n")

rule baseline:
    input:
        conf=dir_work + "baselines/{genome_version}.baseline.samples.list",
        ms=dir_work + "reference/{genome_version}.msisensor.scan.list",
    output:
        dir_work + "baselines/{genome_version}.baseline.tsv"
    run:
        shell(f"{msisensor_pro} baseline  -d {input.ms} -i {input.conf} -o {output} -s 1")
rule run_pro:
    input:
        ms = dir_work + "baselines/{genome_version}.baseline.tsv",
        t=lambda wildcards: samples.loc[wildcards.case, "tumor_path"],
        ref=lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        dir_work + "tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro"
    log:
        dir_work + "tumor_only_output/{case}/{case}.{genome_version}.msisensor-pro.log"
    run:
        shell("{msisensor_pro} pro -d {input.ms} -t {input.t} -g {input.ref}  -o {output} 2>{log}")

rule merge_msi_result:
    input:
        expand(dir_work + "tumor_only_output/{case}/{case}.{{genome_version}}.msisensor-pro",case=samples.index)
    output:
        dir_work + "tumor_only_output.{genome_version}.merge.tsv"
    run:
        output_info = pd.DataFrame(columns=["Total_number_of_sites", "Number_of_unstable_sites", "MSI_score"])
        output_info.index.name="case_name"
        for case in samples.index:
            value_info = [i.split() for i in open(dir_work + f"tumor_only_output/{case}/{case}.{wildcards.genome_version}.msisensor-pro")]
            output_info.loc[case] = value_info[1]
        output_info.to_csv(f"{output}",sep="\t")
