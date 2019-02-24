date_var=`date '+%Y%m%d_%H%M%S'`
dataset_dir=../datasets/$date_var
mkdir $dataset_dir
touch $dataset_dir/baselineTimings.txt

echo "$3" >> $dataset_dir/newEdges.txt

for i in $(seq 0 $3)
do
  printf "\nAdd Edge $i: \n"
  python ../datasets/snap_parse_full.py ../datasets/$1.txt $dataset_dir/$1_$i.metis $dataset_dir/newEdges.edges $2 $i
done

# Overwite firstline of newEdges.edges with number of lines
edge_count="$(wc -l $dataset_dir/newEdges.edges | cut -f1 -d' ')"
echo "${edge_count}"
sed -i "1s/.*/$edge_count/" $dataset_dir/newEdges.edges

mkdir $dataset_dir/baseline
make
for j in $(seq 0 $3)
do
  ./build/Linux-x86_64/kt -kttype=baseline5 -stream=0 $dataset_dir/$1_$j.metis $dataset_dir/baseline/$j.out
  sort $dataset_dir/baseline/$j.out > $dataset_dir/baseline/sorted_$j.out
done

mkdir $dataset_dir/incremental
./build/Linux-x86_64/kt -kttype=baseline5 -stream=1 $dataset_dir/$1_0.metis temp.out
for j in $(seq 1 $3)
do
  sort $dataset_dir/incremental/$j.out > $dataset_dir/incremental/sorted_$j.out
done

awk '{if (NR!=1) n += $1}; END{print "Total Baseline Runtime : " n}' $dataset_dir/baselineTimings.txt
