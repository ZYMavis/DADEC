lreads=sim-ecoli/demo_long.fa
sreads=sim-ecoli/demo_short.fa

../DADEC -s $sreads -l $lreads -o ecoli.fa -S 1 -t 32
