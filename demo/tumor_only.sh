cd tumor_only
../../cpp/msisensor-pro scan -d ../reference/reference.fa -o ../reference/reference.list
../../cpp/msisensor-pro baseline -i baselineBam.configure -d ../reference/reference.list -o baseline
mkdir test_result
../../cpp/msisensor-pro pro -d ./baseline/reference.list_baseline -t test_tumor_sorted.bam -o test_result/test
