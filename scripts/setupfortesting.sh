
cpupower -c all frequency-set --governor performance 

# manually disable turboboost for all cores
NPROCS=`cat /proc/cpuinfo | grep "core id" | wc -l`
NPROCS=$(($NPROCS - 1))
for i in `seq 0 1 $NPROCS`; do
  wrmsr -p $i 0x1a0 0x4000850089;
done

# x86_energy_perf_bias to performance
x86_energy_perf_policy -c all performance
