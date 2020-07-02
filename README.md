[![Bioconda](https://img.shields.io/conda/dn/bioconda/msisensor-pro.svg?label=Bioconda)](https://bioconda.github.io/recipes/msisensor-pro/README.html)
[![Published in Genomics, Proteomics & Bioinformatics](https://img.shields.io/badge/Published%20in-GPB-167DA4.svg)](https://www.sciencedirect.com/science/article/pii/S1672022920300218)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/msisensor-pro/badges/latest_release_date.svg)](https://anaconda.org/bioconda/msisensor-pro)

![MSIsensor-pro](/fig/logo_msisensor-pro.png)

## Please click [here](https://github.com/xjtu-omics/msisensor-pro/wiki) to see more about MSIsensor-pro.

# General introduction

MSIsensor-pro is an updated version of **[msisensor](https://github.com/ding-lab/msisensor)**. MSIsensor-pro evaluates Microsatellite Instability (MSI) for cancer patients with next generation sequencing data. It accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data as input.  MSIsensor-pro introduces a multinomial distribution model to quantify polymerase slippages for each tumor sample and a discriminative sites selection method to enable MSI detection without matched normal samples. For samples of various sequencing depths and tumor purities, MSIsensor-pro significantly outperformed the current leading methods in terms of both accuracy and computational cost. If you want to know more detail about MSIsensor-pro, please see the **[MSIsensor-pro Schematics and Internals](https://github.com/xjtu-omics/msisensor-pro/wiki/MSIsensor-pro-Schematics-and-Internals)** page.

# Scopes of MSIsensor-pro

MSIsensor-pro evaluates MSI status of a given sample with next generation sequencing (NGS) data. If you have normal-tumor paired DNA sequences, you can use _**msi**_ (inherited from msisensor) module to score MSI status while _**pro**_ module would be the option if matched normal is not available.

# How to install MSIsensor-pro?


## [ Install Using Bioconda ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-from-source-code)


      conda install msisensor-pro
      

## [ Install from source code ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-from-source-code)

**( Recommended For Developers )**


#### Install the dependencies
  Dependent packages including zlib, ncurses and nurses-dev are required for MSIsensor-pro. You may already have these prerequisite packages. If not, you need to run the following code to obtain dependent packages.

* For Debian or Ubuntu:

      sudo apt-get install zlib1g-dev libncurses5-dev libncursesw5-dev

* For Fedora, CentOS or RHEL

      sudo yum install zlib-devel ncurses-devel ncurses

#### Build MSIsensor-pro from source code
* colne the repository from our github

      git clone https://github.com/xjtu-omics/msisensor-pro

* make 

      cd msisensor-pro/cpp
      make
 
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


