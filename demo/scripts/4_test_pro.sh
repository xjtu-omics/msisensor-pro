# test msi pro: msi detection with tumor sample
dir_demo=$(realpath ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

echo "================================================================="
echo Test with bam input
${msisensor_pro} pro -d ${dir_demo}/data/data4baseline/baseline_bam/reference.list_baseline \
  -t ${dir_demo}/data/data4test/test_tumor_sorted.bam \
  -o ${dir_demo}/data/output/tumor_only/tumor_only_bam -b 2

echo "================================================================="
echo Test with cram input
${msisensor_pro} pro -d ${dir_demo}/data/data4baseline/baseline_cram/reference.list_baseline \
  -t ${dir_demo}/data/data4test/test_tumor_sorted.cram \
  -g ${dir_demo}/data/reference/reference.fa \
  -o ${dir_demo}/data/output/tumor_only/tumor_only_cram -b 2
