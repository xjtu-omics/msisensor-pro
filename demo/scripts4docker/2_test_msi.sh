# test msi command: msi detection with tumor sample and matched noraml sample
dir_demo=$(realpath  ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}

echo "================================================================="
echo Test with bam input
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro msi \
  -d /data/reference/reference.list \
  -t /data/data4test/test_normal_sorted.bam \
  -n /data/data4test/test_tumor_sorted.bam \
  -o /data/output/tumor_normal/tumor_normal_bam -b 2

echo "================================================================="
echo Test with cram input
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro msi \
  -d /data/reference/reference.list \
  -t /data/data4test/test_normal_sorted.cram \
  -n /data/data4test/test_tumor_sorted.cram \
  -g /data/reference/reference.fa \
  -o /data/output/tumor_normal/tumor_normal_cram -b 2
