# test msi command: msi detection with tumor sample and matched noraml sample
dir_demo=$(realpath  ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

echo "================================================================="
echo Test with bam input
${msisensor_pro} msi -d ${dir_demo}/data/reference/reference.list \
  -t ${dir_demo}/data/data4test/test_normal_sorted.bam \
  -n ${dir_demo}/data/data4test/test_tumor_sorted.bam \
  -o ${dir_demo}/data/output/tumor_normal/tumor_normal_bam -b 2

echo "================================================================="
echo Test with cram input
${msisensor_pro} msi -d ${dir_demo}/data/reference/reference.list \
  -t ${dir_demo}/data/data4test/test_normal_sorted.cram \
  -n ${dir_demo}/data/data4test/test_tumor_sorted.cram \
  -g ${dir_demo}/data/reference/reference.fa \
  -o ${dir_demo}/data/output/tumor_normal/tumor_normal_cram -b 2
