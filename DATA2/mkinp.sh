declare -a Sr

#Sr=( 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.0 )
#Sr=(  0.75 0.80 0.85 0.90 0.95 1.0 )
Sr=( 0.80 0.825 0.85 0.875 0.90 0.925)

declare -i i
echo Sr=${Sr[*]}

N=${#Sr[*]}

i=0
for var in "${Sr[@]}"
do
	echo Sr=${var}
	i=$(($i+1))
	mkdir -p Sr$i
	sed s/DEG_SAT/${var}/ wet0.inp > Sr$i/wet.inp
done
