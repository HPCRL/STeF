

tns=$1
#mode=2

#export OMP_NUM_THREADS=2

for order in 0 1 2 3 4 5 ;
do

	for core in 1 18 ;
	do
		C=$(( $core - 1 ))
		for g in MEM_DP ;
		do
			for mode in 0 1 2 ;
			do
				for saved in 1 2 ;
				do

					cmd="likwid-perfctr -m -g $g -O -C 0-$C ./bin/mttkrp_combined.exe $tns 32 $order $mode $saved"
					echo $cmd
					$cmd
					date
				done
			done
		done
	done

done
