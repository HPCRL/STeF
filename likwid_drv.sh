

tns=~/tensors/nell-2.tns
mode=2

export OMP_NUM_THREADS=2
for g in FLOPS_DP L2 L3 MEM_DP CACHES L2CACHE L3CACHE ;
do
	cmd="likwid-perfctr -m -g $g -O -C 0-27 ./bin/SpTL.exe $tns 32 -1 $mode"
	echo $cmd
	$cmd
done
