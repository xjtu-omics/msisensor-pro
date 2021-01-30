# test scan command: get microsatellite information from reference genome
dir_demo=$(realpath ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}



$msisensor_pro scan -d ${dir_demo}/data/reference/reference.fa -o ${dir_demo}/data/reference/reference.list
