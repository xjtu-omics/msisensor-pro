cd tumor_only
../../msisensor-pro scan -d ../reference/reference.fa -o ../reference/reference.list
../../msisensor-pro baseline -i baselineBam.configure -d ../reference/reference.list -o baseline
../../msisensor-pro pro -d ./baseline/reference.list_baseline -t test_tumor_sorted.bam -o test_result/test
