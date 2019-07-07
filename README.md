# General introduction

MSIsensor-pro is an updated version of **[msisensor](https://github.com/ding-lab/msisensor)**. MSIsensor-pro can evaluate Microsatellite Instability (MSI) for cancer patients with next generation sequencing data and accepts the whole genome sequencing, whole exome sequencing and target region (panel) sequencing data. Given that there are always no normal control samples for many patients who need to do MSI test, MSIsensor-pro adds a new module (pro) for MSI classification with only tumor sequencing data. MSIsensor-pro show a comparable performance with MSIsensor original version. If you want to know more detail about MSIsensor-pro, please see the **[MSIsensor-pro Schematics and Internals](https://github.com/xjtu-omics/msisensor-pro/wiki/MSIsensor-pro-Schematics-and-Internals)** page.

# Scopes of MSIsensor-pro

MSIsensor-pro can evaluate MSI with next generation sequencing (NGS) data. If you have normal-tumor paired DNA sequences, you can use _**msi**_ (inherited from msisensor) module to score MSI status. If you only have tumor sample, you can use _**pro**_ module to do it.

# How to install MSIsensor-pro?




## [ Install from source code ](https://github.com/xjtu-omics/msisensor-pro/wiki/How-to-install-MSIsensor-pro#install-from-source-code)

#### Install the dependencies
   If you want to make MSIsensor-pro by yourself, you need dependent packages including zlib, ncurses and nurses-dev. You may already have these prerequisite packages. If not, you need to run the following code to finish it.

* For Debian or Ubuntu:

     sudo apt-get install zlib1g-dev libncurses5-dev libncursesw5-dev

* For Fedora, CentOS or RHEL

     sudo yum install zlib-devel ncurses-devel ncurses

#### Build MSIsensor-pro from source code
* colne the repository from our github

      git clone https://github.com/xjtu-omics/msisensor-pro

* make 

      cd msisensor-pro
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

  see detail in the **[files format](https://github.com/xjtu-omics/msisensor-pro/wiki/Files-format)** page
## 
