


## _scan_
scan the reference genome to get microsatellites information

       Usage:  msisensor-pro scan [options]

       -d   <string>   reference genome sequences file, *.fasta or *.fa format  [required]
       -o   <string>   output homopolymers and microsatellites file [required] 

       -l   <int>      minimal homopolymer(repeat unit length = 1) size, [default=8]
       -c   <int>      context length (5-32), [default=5]
       -m   <int>      maximal homopolymer size, [default=50]
       -s   <int>      maximal length of microsatellite (1-32) , [default=6]
       -r   <int>      minimal repeat times of microsatellite(repeat unit length>=2), [default=5]
       -p   <int>      output homopolymer only, 0: no; 1: yes, [default=0]
       -h   help

## msi
 evaluate MSI using paired tumor-normal sequencing data (msisensor)

       Usage:  msisensor-pro msi [options]

       -d   <string>   homopolymers and microsatellites file [required]
       -n   <string>   normal bam file with index [required]
       -t   <string>   tumor  bam file with index [required]
       -g   <string>   reference file [required if *.cram for -t]
       -o   <string>   output path (Ending with a slash is not allowed.) [required] 

       -f   <double>   FDR threshold for somatic sites detection, [default=0.05]
       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, [default=15]
       -z   <int>      coverage normalization for paired tumor and normal data, 0: no; 1: yes, [default=0]
       -p   <int>      minimal homopolymer size for distribution analysis, [default=8]
       -m   <int>      maximal homopolymer size for distribution analysis, [default=50]
       -s   <int>      minimal microsatellite size for distribution analysis, [default=5]
       -w   <int>      maximal microsatellite size for distribution analysis, [default=40]
       -u   <int>      span size around window for extracting reads, [default=500]
       -b   <int>      threads number for parallel computing, [default=1]
       -x   <int>      output homopolymer only, 0: no; 1: yes, [default=0]
       -y   <int>      output microsatellites only, 0: no; 1: yes, [default=0]
       -0   <int>      output site have no read coverage, 1: no; 0: yes, [default=1]
       -h   help


## _pro_
   evaluate MSI using single (tumor) sample sequencing data

       Usage:  msisensor-pro pro [options]

       -d   <string>   homopolymers and microsatellites file [required] 
       -t   <string>   bam/cram file of tumor/normal(for baseline building) sample [required] 
       -g   <string>   reference file [required if *.cram for -t]
       -o   <string>   output path (Ending with a slash is not allowed.) [required]

       -i   <double>   minimal threshold for instable sites detection (just for tumor only data), [default=0.1]
       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, [default=15]
       -p   <int>      minimal homopolymer size for distribution analysis, [default=8]
       -m   <int>      maximal homopolymer size for distribution analysis, [default=50]
       -s   <int>      minimal microsatellite size for distribution analysis, [default=5]
       -w   <int>      maximal microsatellite size for distribution analysis, [default=40]
       -u   <int>      span size around window for extracting reads, [default=500]
       -b   <int>      threads number for parallel computing, [default=1]
       -x   <int>      output homopolymer only, 0: no; 1: yes, [default=0]
       -y   <int>      output microsatellite only, 0: no; 1: yes, [default=0[
       -0   <int>      output site have no read coverage, 1: no; 0: yes, [default=1]
       -h   help

## _baseline_
build baseline for tumor only detection


       Usage:  msisensor-pro baseline [options]

       -d   <string>   homopolymer and microsatellite file [required]
       -i   <string>   configure files for building baseline (text file) [required]
            you need to provide the output (*_all) from pro command 
            e.g.
             ----------------------------------
              case1	/path/to/case1_sorted_all 
              case2	/path/to/case2_sorted_all 
              case3	/path/to/case3-sorted_all 
             ----------------------------------
       -o   <string>   output path for baseline [required] 

       -s   <double>   microsatellite sites with support from fewer than -d samples will not pass quality control, [default=10]
       -h   help