## Best Practices for MSI evaluation using MSIsensor-pro
We recommend using Snakemake to run the MSIsensor-pro pipeline,
especially when you have a large sample size. For both tumor-only and tumor-normal modes, 
please first modify the parameters for running MSIsensor-pro in the config.yaml and *.smk files. 
Then, run your pipeline using Snakemake.




# 1. configure the config.yaml file

 * ## config.yaml of tumor-normal paired model.

    ```yaml
    msisensor-pro: # path of the msisensor-pro
      /path/of/msisensor-pro/cpp/msisensor-pro
    dir_work: # path of your result output
      /path/of/your/pipeline_output/tumor_normal/
    reference: # path of your reference genome
      /path/of/reference/reference.fa
    genome_version: # version ID 
      GRCh38
    sample_path: #path of the sample path, see the  template at demo/data_input/data4test/tumor_normal_paired.samples.info.tsv
      /path/of/tumor_normal_paired/tumor_normal_paired.samples.info.tsv
    
    ```

* ## config.yaml of tumor-only paired model.
    ```yaml
     msisensor-pro: # path of the msisensor-pro
       /path/of/msisensor-pro/cpp/msisensor-pro
     dir_work:
      /path/of/your/pipeline_output/tumor_only/
     reference: # path of your reference genome
      /path/of/reference/reference.fa
     genome_version:
       GRCh38
     samples4baseline: #path of the sample info for baseline building, see template at demo/data_input/data4baseline/sample4baseline_info.tsv
       /path/of/tumor_only/data4baseline/sample4baseline_info.tsv
     samples4test: #path of the sample info for MSI evaluating, see template at demo/data_input/data4test/tumor_only.samples.info.tsv
       /path/of/tumor_only/data4test/tumor_only.samples.info.ts
     ```
 
