# test scan: get microsatellite infomation from reference
dir_demo=`readlink  -f ..`
msisensor_pro=`readlink -f ../../binary/msisensor-pro`
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}
${msisensor_pro} scan -d ${dir_demo}/data/reference/reference.fa -o ${dir_demo}/data/reference/reference.list
