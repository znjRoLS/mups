if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Using default number of proc"
    mpirun bin/dz2z1 1000 0 1000
    mpirun bin/dz2z1 1000 100 10000
    mpirun bin/dz2z1 1000000 0 1000
    mpirun bin/dz2z1 1000000 1000 10000
    mpirun bin/dz2z1 1000000000 0 1000
    mpirun bin/dz2z1 1000000000 100 10000
 else
    echo "Running with $1 number of procesors"
    mpirun -np $1 bin/dz2z1 1000 0 1000
    mpirun -np $1 bin/dz2z1 1000 100 10000
    mpirun -np $1 bin/dz2z1 1000000 0 1000
    mpirun -np $1 bin/dz2z1 1000000 1000 10000
    mpirun -np $1 bin/dz2z1 1000000000 0 1000
    mpirun -np $1 bin/dz2z1 1000000000 100 10000
fi

