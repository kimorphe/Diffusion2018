declare -a Sr
declare -i seed

seed=-3
Sr=( 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
declare -i i
#echo Sr=${Sr[*]}

N=${#Sr[*]}

i=0
for var in "${Sr[@]}"
do
	echo Sr=${var}
	i=$(($i+1))
	mkdir -p Sr$i
	sed s/DEG_SAT/${var}/ wet0.inp > tmp.inp
	sed s/SEED/${seed}/ tmp.inp > Sr$i/wet.inp
done
rm tmp.inp
