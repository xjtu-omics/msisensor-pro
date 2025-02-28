[![Published in Genomics, Proteomics & Bioinformatics](https://img.shields.io/badge/Published%20in-GPB-167DA4.svg)](https://www.sciencedirect.com/science/article/pii/S1672022920300218)
![GitHub last commit](https://img.shields.io/github/last-commit/xjtu-omics/msisensor-pro)
[![GitHub Release Date](https://img.shields.io/github/release-date/xjtu-omics/msisensor-pro)](https://github.com/xjtu-omics/msisensor-pro/releases)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/xjtu-omics/msisensor-pro?include_prereleases)](https://github.com/xjtu-omics/msisensor-pro/releases)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/pengjia1110/msisensor-pro)](https://hub.docker.com/repository/docker/pengjia1110/msisensor-pro)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/msisensor-pro.svg?label=Bioconda)](https://bioconda.github.io/recipes/msisensor-pro/README.html)
![Docker Pulls](https://img.shields.io/docker/pulls/pengjia1110/msisensor-pro)
![GitHub all releases](https://img.shields.io/github/downloads/xjtu-omics/msisensor-pro/total?label="Github") 
![GitHub Repo stars](https://img.shields.io/github/stars/xjtu-omics/msisensor-pro)


![MSIsensor-pro](./fig/logo_msisensor-pro.png)




[//]: # (## Please click [here]&#40;https://github.com/xjtu-omics/msisensor-pro/wiki&#41; to see more about MSIsensor-pro in Wiki. )

## Welcome to try our new software [MSIsensor-RNA](https://github.com/xjtu-omics/msisensor-rna) for MSI detection with RNA-seq data!


# Introduction of MSIsensor-pro
MSIsensor-pro is an updated version of **[msisensor](https://github.com/ding-lab/msisensor)**.
MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with next generation sequencing
data. It accepts the whole genome sequencing, whole exome sequencing and target region (panel)
sequencing data as input.  MSIsensor-pro introduces a multinomial distribution model
to quantify polymerase slippages for each tumor sample and a discriminative sites selection
method to enable MSI detection without matched normal samples. For samples of various
sequencing depths and tumor purities, MSIsensor-pro significantly outperformed
the current leading methods which required matched normal samples in terms of both accuracy
and computational cost. If you want to know more detail about MSIsensor-pro, please see
our research paper at [_**Genomics, Proteomics & Bioinformatics**_](https://www.sciencedirect.com/science/article/pii/S1672022920300218) page.

## [Best practices pipeline](./snakemake)
## [How to install MSIsensor-pro?](./docs/3_Installation.md)
## [Main commands of MSIsensor-pro](./docs/4_command.md)
## [License of MSIsensor-pro](./docs/2_License.md)
## [File format of MSIsensor-pro](./docs/5_files_type.md)
## [Frequently asked questions](./docs/6_Frequently_asked_questions.md)



# Contact

* If you have any technical questions, please [open an issue](https://github.com/xjtu-omics/msisensor-pro/issues/new/choose). If you don't get a prompt response(maybe two working day), please contact with Peng Jia (pengjia@xjtu.edu.cn).

* If you want to apply the MSIsensor-pro to commercial purposes, please click [here](./docs/2_License.md).



# Citation
Peng Jia, Xiaofei Yang, Li Guo, Bowen Liu, Jiadong Lin, Hao Liang, et al. 
MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability. 
Genomics Proteomics Bioinformatics 2020,18(1).  [PDF](https://www.sciencedirect.com/science/article/pii/S1672022920300218)



See more details in the **[Key Commands](https://github.com/xjtu-omics/msisensor-pro/wiki/Key-Commands)** page and **[Best Practices](https://github.com/xjtu-omics/msisensor-pro/wiki/Best-Practices)** page.




# Contributors
* Peng Jia 
* Bowen Liu 
* Hao Liang 
* Mingzhe Duan


# Disclaimer
MSIsensor-pro, its developers, and owners do not intend to provide medical advice nor to replace the medical advice or treatment from your doctor. This information is not intended to diagnose, treat, cure, or prevent any disease.
