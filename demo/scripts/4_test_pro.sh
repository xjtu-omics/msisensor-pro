# test msi pro: msi detection with tumor sample
dir_demo=$(readlink -f ..)
msisensor_pro=$(readlink -f ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

${msisensor_pro} pro -d ${dir_demo}/data/data4baseline/baseline/reference.list_baseline \
  -t ${dir_demo}/data/data4test/test_normal_sorted.bam \
  -o ${dir_demo}/data/output/tumor_only/tumor_only
