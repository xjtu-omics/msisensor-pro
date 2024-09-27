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


dir_work = config["dir_work"]
dir_work = dir_work if dir_work.endswith("/") else dir_work + "/"
samples = pd.read_table(f"{config['sample_path']}",index_col=0)

msisensor_pro = config["msisensor-pro"]
genome_version = config["genome_version"]

genome_reference = {config["genome_version"]: config["reference"]}


rule all:
    input:
        dir_work + f"tumor_normal_output.{genome_version}.merge.tsv"


rule scan:
    input:
        lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        dir_work + "reference/{genome_version}.msisensor.scan.list"
    log:
        dir_work + "reference/{genome_version}.msisensor.scan.log"
    run:
        shell("{msisensor_pro} scan -d {input} -o {output} 2>{log} 1>{log}")

rule msisensor_msi:
    input:
        t=lambda wildcards: samples.loc[wildcards.case, "tumor_path"],
        n=lambda wildcards: samples.loc[wildcards.case, "normal_path"],
        ms=dir_work + "reference/{genome_version}.msisensor.scan.list",
        ref=lambda wildcards: genome_reference[wildcards.genome_version]
    output:
        dir_work + "tumor_normal_output/{case}/{case}.{genome_version}.msisensor"
    log:
        dir_work + "tumor_normal_output/{case}/{case}.{genome_version}.msisensor.log"
    run:
        shell("{msisensor_pro} msi -d {input.ms} -n {input.n} -t {input.t} -g {input.ref}  -o {output} 2>{log}")

rule merge_msi_result:
    input:
        expand(dir_work + "tumor_normal_output/{case}/{case}.{{genome_version}}.msisensor",case=samples.index)
    output:
        dir_work + "tumor_normal_output.{genome_version}.merge.tsv"
    run:
        output_info = pd.DataFrame(columns=["Total_number_of_sites", "Number_of_unstable_sites", "MSI_score"])
        output_info.index.name="case_name"
        for case in samples.index:
            value_info = [i.split() for i in open(dir_work + f"tumor_normal_output/{case}/{case}.{wildcards.genome_version}.msisensor")]
            output_info.loc[case] = value_info[1]
        output_info.to_csv(f"{output}",sep="\t")
# print(value_info)
