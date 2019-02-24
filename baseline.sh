# Parse input arguments
dataset_name=$1
initial_num_edges=$2
if [ -z "$3" ]; then
  streaming_edges="all"
else
  streaming_edges=$3
fi
dataset_dir=../datasets/$1_$2_$streaming_edges

# Create essential files
mkdir $dataset_dir
mkdir $dataset_dir/inputs

# Create all .metis files
printf "\nCreating baseline files..."
python ../datasets/snap_parse_baseline.py ../datasets/$dataset_name $dataset_dir/inputs/ $dataset_dir/inputs/newEdges $initial_num_edges $streaming_edges

# Overwite firstline of newEdges.edges with number of lines
edge_count="$(wc -l $dataset_dir/inputs/newEdges.edges | cut -f1 -d' ')"
echo "${edge_count}"
sed -i "1s/.*/$edge_count/" $dataset_dir/inputs/newEdges.edges

# Run k-truss for all datasets individually
make
mkdir $dataset_dir/outputs_baseline
for i in $(seq 0 $edge_count)
do
  ./build/Linux-x86_64/kt -kttype=baseline5 -stream=0 $dataset_dir/inputs/$i.metis $dataset_dir/outputs_baseline/$i.out
done
