# Parse input arguments
date_var=`date '+%Y%m%d_%H%M%S'`
dataset_directory=$1
dataset_dir=../datasets/$dataset_directory

# Make necessary directories/files
output_dir=$dataset_dir/incremental_$date_var
mkdir $output_dir
cp $dataset_dir/inputs/0.metis $output_dir
cp $dataset_dir/inputs/newEdges.edges $output_dir

# Stream
make
# valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./build/Linux-x86_64/kt -kttype=baseline5 -stream=1 $output_dir/0.metis temp.out
./build/Linux-x86_64/kt -kttype=baseline5 -stream=1 $output_dir/0.metis temp.out

# Compare Baseline and Incremental timings
awk '{if (NR!=1) n += $1}; END{print "Total Baseline Runtime : " n}' $dataset_dir/inputs/baselineTimings.txt
awk '{if (NR!=1) n += $1}; END{print "Total Baseline Ops : " n}' $dataset_dir/inputs/baselineOps.txt
# awk '{n += $1}; END{print "Total Incremental Runtime : " n}' $output_dir/incrementalTimings.txt
