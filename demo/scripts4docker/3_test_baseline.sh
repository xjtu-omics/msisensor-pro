# test baseline command: build baseline for tumor only methods
dir_demo=$(realpath ..)
msisensor_pro=$(realpath ../../binary/msisensor-pro)
echo demo_dir: ${dir_demo}
echo msisensor-pro: ${msisensor_pro}

echo "================================================================="
echo Test with bam input
#cd ${dir_demo}/data/data4baseline
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest  msisensor-pro baseline \
  -i /data/data4baseline/baseline.bam.configure4docker \
  -d /data/reference/reference.list \
  -o /data/data4baseline/baseline_bam -b 2 


echo "================================================================="
echo Test with cram input
docker run -v ${dir_demo}/data:/data  pengjia1110/msisensor-pro:latest msisensor-pro baseline \
  -i /data/data4baseline/baseline.cram.configure4docker \
  -d /data/reference/reference.list \
  -g /data/reference/reference.fa \
  -o /data/data4baseline/baseline_cram -b 2
