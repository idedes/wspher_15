# My job

source /libcern/intel/Compiler/2011.9.293/setup.sh 

tar zxvf wspher15.tar.gz

./wspher_15_split
./wspher_15_cmp_ifc

./wspher15_exe < wspher_15.d > wspher_15.out

tar zcvf output.tar.gz DensitCurves EnergyLevels LogFiles wspher_15.out

