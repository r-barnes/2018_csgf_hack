sshfs -o ServerAliveInterval=15,ServerAliveCountMax=3 -C spinyfan@home.ccs.ornl.gov:/ccs/home/spinyfan ~/temp/ornl

module switch PrgEnv-pgi PrgEnv-gnu
module load scorep
CXX=scorep-CC
CXXFLAGS=-O3 -g
make


qsub -I -X -A CSC261 -q debug -l nodes=1,walltime=30:00
cd $MEMBERWORK/csc261
cp ~/path/to/gppKerSeq.ex ./