for i in $(seq 1 $2)
do
  sort $1/../outputs_baseline/$i.out > $1/../outputs_baseline/sorted$i.out
  sort $1/$i.out > $1/sorted$i.out
  diff $1/../outputs_baseline/sorted$i.out $1/sorted$i.out
done
