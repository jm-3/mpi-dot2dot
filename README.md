# mpi-dot2dot

MPI implementation of the algorithm Dot2Dot which is a method for tandem repeats discovery within DNA sequences

## Contents:

- datasets: Due to its size, datasets are not included. Instead, a script is provided to download them
- doc: Dot2Dot user manual
- scripts: useful scripts to help to automatize the experiments
    - slurm: scripts related to this job manager
        - slanzarv: creates the submit script
        - benchmarks.sh: sends jobs to cluster, using slanzarv. Results are stored in path: preffixDir/datasetName/algorithm_version/num_procs/test_number/ . 3 test are performed per each configuration. 
        - getSacctValues.sh: extracts job metrics from sacct. Automatically, it searches recursively for jobs ID and extract the pertinent info
		- plot: plotting results python scripts
- src: source code
	- Dot-1.0.p3: initial Dot2Dot code. Available from https://github.com/Gege7177/Dot2dot
	- secuencial: Dot2Dot sequential version. Threads removed
	- mpi: multiprocessor version (MPI). Block and balanced Scheduler
	- hybrid: hybrid version (MPI + OpenMP)
	- hybrid-variable: hybrid version with dynamic resource allocation		
- tests: benchmarks and other tests about performance
	- lustre: stripe size and count testing over big files
	- secuencias: counting multisequence files tests
	- pthreads: alternative pthreads implentation
