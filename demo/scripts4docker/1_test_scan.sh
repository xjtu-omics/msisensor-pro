# test scan command: get microsatellite information from reference genome
dir_demo=$(realpath ..)
echo demo_dir: ${dir_demo}



docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro  scan -d /data/reference/reference.fa -o data/reference/reference.list
