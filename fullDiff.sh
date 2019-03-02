for i in $(seq 1 $2)
do
  diff $1/baseline/sorted_$i.out $1/sorted_$i.out
done
