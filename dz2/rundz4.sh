if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Using default number of proc"
    mpirun bin/dz2z4 200 500 0.1 sol4_200x500.txt
    mpirun bin/dz2z4 500 500 0.1 sol4_500x500.txt
 else
    echo "Running with $1 number of procesors"
    mpirun -np $1 bin/dz2z4 200 500 0.1 sol4_200x500.txt
    mpirun -np $1 bin/dz2z4 500 500 0.1 sol4_500x500.txt
fi