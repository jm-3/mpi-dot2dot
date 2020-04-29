# mpi-dot2dot

MPI implementation of the algorithm Dot2Dot which is a method for tandem repeats discovery within DNA sequences

## Contents:

- Dot-1.0.p3: initial Dot2Dot code. Available from https://github.com/Gege7177/Dot2dot 
- Seq-Dot2Dot: sequential version of the original Dot2Dot, obtained from removing parallel parts (POSIX threads and lockers)
- MPI-Dot2Dot: MPI version of Dot2Dot, developed from the sequential one
- doc: Dot2Dot user manual