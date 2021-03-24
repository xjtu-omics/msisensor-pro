# test msi pro: msi detection with tumor sample
dir_demo=$(realpath ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}

echo "================================================================="
echo Test with bam input
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro pro \
  -d /data/data4baseline/baseline_bam/reference.list_baseline \
  -t /data/data4test/test_tumor_sorted.bam \
  -o /data/output/tumor_only/tumor_only_bam -b 2

echo "================================================================="
echo Test with cram input
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro pro \
  -d /data/data4baseline/baseline_cram/reference.list_baseline \
  -t /data/data4test/test_tumor_sorted.cram \
  -g /data/reference/reference.fa \
  -o /data/output/tumor_only/tumor_only_cram -b 2
