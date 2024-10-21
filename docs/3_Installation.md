# How to install MSIsensor-pro?

## Directly using binary version 
            
      wget https://github.com/xjtu-omics/msisensor-pro/raw/master/binary/msisensor-pro[v1.x.x]
      chmod +x msisensor-pro[v1.x.x]
      export PATH=`pwd`:$PATH

If some dynamic libraries are missing, you need to install them and add their paths to the LD_LIBRARY_PATH variable.

## Install Using Docker 

       docker pull pengjia1110/msisensor-pro   
       docker run pengjia1110/msisensor-pro msisensor-pro

When running msisensor-pro with Docker, please note that you need to map the local path to the internal Docker path using the -v option.

## Install Using Bioconda


      conda install msisensor-pro=v1.2.0 # install v1.2.0 version of msisensor-pro with bioconda 
      conda install msisensor-pro # install latest version of msisensor-pro with bioconda 



## Install from source code

**( Recommended For Developers )**


#### Install the dependencies
Dependent packages including zlib, ncurses and nurses-dev are required for MSIsensor-pro. You may already have these prerequisite packages. If not, you need to run the following code to obtain dependent packages.

* For Debian or Ubuntu:

      sudo apt-get install libbz2-dev zlib1g-dev libcurl4-openssl-dev libncurses5-dev libncursesw5-dev

* For Fedora, CentOS or RHEL

      sudo yum install bzip2-devel xz-devel zlib-devel ncurses-devel ncurses

#### Build MSIsensor-pro from source code
* clone the repository from our github

      git clone https://github.com/xjtu-omics/msisensor-pro

* make

      cd msisensor-pro/
      ./INSTALL

* install

      sudo mv binary/msisensor-pro /usr/local/bin/

  You can also choose to add the path of msisensor-pro to the PATH variable or run it directly using the absolute path.

