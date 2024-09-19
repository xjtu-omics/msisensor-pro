# ======================================================================================================================
# Project: Project_ABC
# Script : filter_fann.smk TODO check
# Author : Peng Jia
# Date   :  2024/5/27
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
configfile: "baseline.conf.yaml"
import pandas as pd

print(config)
dir_work = config["dir_work"]
dir_work = dir_work if dir_work.endswith("/") else dir_work + "/"
baseline_bams = config["samples4baseline"]
samples = pd.read_table(f"{baseline_bams}",header=None,)
samples.columns = ["sample", "path"]
samples.index = samples["sample"]
msisensor_pro = config["msisensor-pro"]
genome_version = config["genome_version"]
# print(samples)
print(dir_work + f"baselines/{genome_version}.baseline.file")

rule all:
    input:
        dir_work + f"baselines/{genome_version}.baseline.file"
        # dir_work + f"baselines/{genome_version}.baseline.file"

rule pro_pre:
    input:
        bam=lambda wildcards: samples.loc[wildcards.sample, "path"],
        ms=config["ms_list"],
        ref=config["reference"],
    output:
        dir_work + "baselines/details/{sample}.{genome_version}.baseline.out"
    threads: 2
    run:
        shell("{msisensor_pro} pro -d {input.ms} -t {input.bam} -g {input.ref} -o {output}")

rule merge:
    input:
        expand(dir_work + "baselines/details/{sample}.{{genome_version}}.baseline.out",sample=samples["sample"])
    output:
        dir_work + "baselines/{genome_version}.baseline.samples.list"
    run:
        with open(f"{output}","w") as f:
            for i in input:
                f.write(f"{i.split('/')[-1]} \t{i}_all\n")


rule baseline:
    input:
        conf=dir_work + "baselines/{genome_version}.baseline.samples.list",
        ms=config["ms_list"],

    output:
        dir_work + "baselines/{genome_version}.baseline.file"
    run:
        shell(f"{msisensor_pro} baseline  -d {input.ms} -i {input.conf} -o {output} -s 1")

