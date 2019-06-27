#pore
Sr=( 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
sd=( 0 -1 -2 -3 -4 -5 -6 -7 -8 -9)
n=${#Sr[*]}

i=0
for var in "${Sr[@]}" 
do
	i=$(($i+1))
	cd Sr$i
	echo $PWD 
	j=0
	for seed in "${sd[@]}"
	do 
		j=$(($j+1))

		sed s/DEG_SAT/${var}/ ../wet0.inp > tmp.inp
		sed s/SEED/${seed}/ tmp.inp > wet.inp

		../../wet
		../../rwk ../rwk.inp
		if ((${i}==9)) 
		then
			exit 0	
		fi

		mkdir -p SEED$j
		mv rwk*.out SEED$j
		mv pore.dat SEED$j
		mv u2b.out SEED$j
		mv wet.inp SEED$j
	done
	cd ../
done

