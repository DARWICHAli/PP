export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/a.txt out/a_best.out
export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/b.txt out/b_best.out
#export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/blong.txt out/blong_best.out
export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/c.txt out/c_best.out
#export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/clong.txt out/clong_best.out
export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/d.txt out/d_best.out
#export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/dlong.txt out/dlong_best.out
export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/e.txt out/e_best.out
export OMP_NUM_THREADS=$1; mpirun -n $2 ./checker input/f.txt out/f_best.out
