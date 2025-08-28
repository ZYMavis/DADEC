lreads=/home/yczhang/zyc/final_result/ecoil/data/long/long_reads.fa
sreads=/home/yczhang/zyc/final_result/ecoil/data/short/short_reads.fa

DADEC -s $sreads -l $lreads -o ecoli.fa -S 1 -t 32
