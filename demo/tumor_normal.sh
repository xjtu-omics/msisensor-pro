cd tumor_normal
../../cpp/msisensor-pro scan -d ../reference/reference.fa -o ../reference/reference.list
../../cpp/msisensor-pro msi -d ../reference/reference.list -t test_tumor_sorted.bam -n test_normal_sorted.bam -o test
