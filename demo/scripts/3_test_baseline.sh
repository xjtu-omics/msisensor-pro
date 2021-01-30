# test baseline command: build baseline for tumor only methods
dir_demo=$(realpath ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

echo "================================================================="
echo Test with bam input
cd ${dir_demo}/data/data4baseline
${msisensor_pro} baseline -i baseline.bam.configure \
  -d ${dir_demo}/data/reference/reference.list \
  -o baseline_bam -b 2
cd -

echo "================================================================="
echo Test with cram input
cd ${dir_demo}/data/data4baseline
${msisensor_pro} baseline -i baseline.cram.configure \
  -d ${dir_demo}/data/reference/reference.list \
  -g ${dir_demo}/data/reference/reference.fa \
  -o baseline_cram -b 2
cd -
