
[![Published in Genomics, Proteomics & Bioinformatics](https://img.shields.io/badge/Published%20in-GPB-167DA4.svg)](https://www.sciencedirect.com/science/article/pii/S1672022920300218)
![GitHub last commit](https://img.shields.io/github/last-commit/xjtu-omics/msisensor-pro)
[![GitHub Release Date](https://img.shields.io/github/release-date/xjtu-omics/msisensor-pro)](https://github.com/xjtu-omics/msisensor-pro/releases)
[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/xjtu-omics/msisensor-pro?include_prereleases)](https://github.com/xjtu-omics/msisensor-pro/releases)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/pengjia1110/msisensor-pro)](https://hub.docker.com/repository/docker/pengjia1110/msisensor-pro)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/msisensor-pro.svg?label=Bioconda)](https://bioconda.github.io/recipes/msisensor-pro/README.html)
![Docker Pulls](https://img.shields.io/docker/pulls/pengjia1110/msisensor-pro)
![GitHub all releases](https://img.shields.io/github/downloads/xjtu-omics/msisensor-pro/total?label="Github")


![MSIsensor-pro](https://raw.githubusercontent.com/xjtu-omics/msisensor-pro/master/fig/logo_msisensor-pro.png)


## Please click [here](https://github.com/xjtu-omics/msisensor-pro/wiki) to see more about MSIsensor-pro in Wiki. 

## Welcome to try our new software [MSIsensor-RNA](https://github.com/xjtu-omics/msisensor-rna) for MSI detection with RNA-seq data!


# Contact

If you want to apply the MSIsensor-pro to commercial purposes, 
please contact Peng Jia (pengjia@stu.xjtu.edu.cn) or 
Kai Ye (kaiye@xjtu.edu.cn) for a license
and get more services.

# License

MSIsensor-pro is free for non-commercial use
by academic, government, and non-profit/not-for-profit institutions. A
commercial version of the software is available and licensed through
Xi’an Jiaotong University. For more information, please contact with
Peng Jia (pengjia@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).



# Citation
Peng Jia, Xiaofei Yang, Li Guo, Bowen Liu, Jiadong Lin, Hao Liang, et al. MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability. Genomics Proteomics Bioinformatics 2020,18(1).  [PDF](https://www.sciencedirect.com/science/article/pii/S1672022920300218)


# General introduction

MSIsensor-pro is an updated version of **[msisensor](https://github.com/ding-lab/msisensor)**. MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with next generation sequencing data. It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input.  MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample and a discriminative sites selection method to enable MSI detection without matched normal samples. For samples of various sequencing depths and tumor purities, MSIsensor-pro significantly outperformed the current leading methods in terms of both accuracy and computational cost. If you want to know more detail about MSIsensor-pro, please see the **[MSIsensor-pro Schematics and Internals](https://github.com/xjtu-omics/msisensor-pro/wiki/MSIsensor-pro-Schematics-and-Internals)** page.

# Scopes of MSIsensor-pro

MSIsensor-pro evaluates MSI status of a given sample with next generation sequencing (NGS) data. If you have normal-tumor paired DNA sequences, you can use _**msi**_ (inherited from msisensor) module to score MSI status while _**pro**_ module would be the option if matched normal is not available.

# How to install MSIsensor-pro?

## [ Directly using binary version ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#directly-using-binary-version) 

      wget https://github.com/xjtu-omics/msisensor-pro/raw/master/binary/msisensor-pro
      chmod +x msisensor-pro 
      export PATH=`pwd`:$PATH

## [ Install Using Docker ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-using-docker)

       docker pull pengjia1110/msisensor-pro   
       docker run pengjia1110/msisensor-pro msisensor-pro

## [ Install Using Bioconda ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-from-source-code)


      conda install msisensor-pro
      

## [ Install from source code ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-from-source-code)

**( Recommended For Developers )**


#### Install the dependencies
  Dependent packages including zlib, ncurses and nurses-dev are required for MSIsensor-pro. You may already have these prerequisite packages. If not, you need to run the following code to obtain dependent packages.

* For Debian or Ubuntu:

      sudo apt-get install libbz2-dev zlib1g-dev libcurl4-openssl-dev libncurses5-dev libncursesw5-dev

* For Fedora, CentOS or RHEL

      sudo yum install bzip2-devel xz-devel zlib-devel ncurses-devel ncurses

#### Build MSIsensor-pro from source code
* colne the repository from our github

      git clone https://github.com/xjtu-omics/msisensor-pro

* make 

      cd msisensor-pro/
      ./INSTALL
 
* install

      sudo mv msisensor-pro /usr/local/bin/


 

# How to use MSI ? 

### Usage:   
   
      msisensor-pro <command> [options]

### Key Commands:

* **scan**
	  
      scan the reference genome to get microsatellites information

* **baseline**

	   build baseline for tumor only detection

* **msi**

	   evaluate MSI using paired tumor-normal sequencing data

* **pro**

	   evaluate MSI using single (tumor) sample sequencing data 

See more detail in the **[Key Commands](https://github.com/xjtu-omics/msisensor-pro/wiki/Key-Commands)** page and **[Best Practices](https://github.com/xjtu-omics/msisensor-pro/wiki/Best-Practices)** page.

# Files  format

  see details in the **[Files format](https://github.com/xjtu-omics/msisensor-pro/wiki/Files-format)** page in WiKi.
## 
# Frequently asked questions
   
  see details in the **[Frequently asked questions](https://github.com/xjtu-omics/msisensor-pro/wiki/Frequently-Asked-Questions)** page in WiKi.


# Citation
  Peng Jia, Xiaofei Yang, Li Guo, Bowen Liu, Jiadong Lin, Hao Liang, et al. MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability. Genomics Proteomics Bioinformatics 2020,18(1).  [PDF](https://www.sciencedirect.com/science/article/pii/S1672022920300218)   
   
# Contact

If you have any questions, please contact with Peng Jia (pengjia@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).

# Contributors
* Peng Jia 
* Bowen Liu 
* Hao Liang 
* Mingzhe Duan

