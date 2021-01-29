# test msi command: msi detection with tumor sample and matched noraml sample
dir_demo=$(readlink -f ..)
msisensor_pro=$(readlink -f ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

${msisensor_pro} msi -d ${dir_demo}/data/reference/reference.list \
  -t ${dir_demo}/data/data4test/test_normal_sorted.bam \
  -n ${dir_demo}/data/data4test/test_tumor_sorted.bam \
  -o ${dir_demo}/data/output/tumor_normal/tumor_normal
