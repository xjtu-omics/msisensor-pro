# clean all generated data
dir_demo=$(readlink -f ..)
rm -f ${dir_demo}/data/reference/reference.list
rm -rf ${dir_demo}/data/output/tumor_normal/*
rm -rf ${dir_demo}/data/data4baseline/baseline
rm -rf ${dir_demo}/data/output/tumor_only/*
echo "clean all generated data"
