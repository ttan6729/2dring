## 2dring
#2dring - What is it?
2dring is a utility designed for all-pairs correlation coefficient calcualtion in distributed computing (also potential for other all-pairs calculation). It reduces a half communcation cost compares to traditional all-pairs calcualtion scheme and hence sutiable for computing in HPC cluster.


Usage: mpirun -np <numprocs> 2dring <N> [number of vectors] <M> [length of each vector] <output> [path of output] <src> [path of src] <x_d> [x dimension of 2dring]
example: mpirun -np 32 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2 or srun -n 16 -d 24 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2
Notice, currently numprocs should be divisible by xd, and numprocs/xd should be even\n
