# test baseline command: build baseline for tumor only methods
dir_demo=$(readlink -f ..)
msisensor_pro=$(readlink -f ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

cd ${dir_demo}/data/data4baseline
${msisensor_pro} baseline -i baselineBam.configure \
  -d ${dir_demo}/data/reference/reference.list \
  -o baseline
cd -
