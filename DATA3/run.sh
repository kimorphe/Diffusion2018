#pore
Sr=( 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
n=${#Sr[*]}

i=0
for var in "${Sr[@]}" 
do
	i=$(($i+1))
	cd Sr$i
	echo $PWD 
	../../wet 
	../../rwk ../rwk.inp
	cd ../
done

