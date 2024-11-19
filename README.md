# MPI Mandelbrot Set Calculation (Manager-Worker Pattern)

This program calculates the Mandelbrot set using a parallelized manager-worker approach with MPI. The computation is distributed across multiple processes, where the manager delegates tasks to worker processes, which compute and return results incrementally.

## Getting Started

### Prerequisites

- **MPI Library**: Ensure MPI (such as OpenMPI) is installed on your system.
- **C Compiler**: You'll need a C compiler, such as `gcc`.

### Compilation

To compile the program, run:

```bash
mpicc -o mandelbrot_mw mandelbrot_mw.c -lm
```

### Running the Program

To execute the program, use the provided SLURM script or manually run it with MPI:

```bash
sbatch run_mandelbrot_mw.sh
```

Alternatively:

```bash
mpirun -np <num_processes> ./mandelbrot_mpi_mw
```

Replace `<num_processes>` with the number of parallel processes you wish to use.

## Program Details

### Overview

This program divides the complex plane into sections, assigning rows of calculations to worker processes. Each process computes iterations for the Mandelbrot equation and sends results back to the manager for aggregation.

- **Manager**: Coordinates tasks and collects results from workers.
- **Workers**: Compute rows of the Mandelbrot set and send results to the manager incrementally.

### Key Components

- **calc_vals**: Computes the number of iterations for points in the Mandelbrot set.
- **do_communication**: Originally used for final data collation (now replaced with incremental communication).
- **write_to_file**: Saves the results to a file (`mandelbrot.dat`).

### Output

The results are saved in `mandelbrot.dat` (real and imaginary coordinates with corresponding iteration counts) if `doIO` is enabled. Performance stats (execution time and process count) are printed to the console by the manager.

## Customization

Modify the following constants to adjust the calculation's precision or domain:

- `N_RE`: Number of intervals on the real axis.
- `N_IM`: Number of intervals on the imaginary axis.
- `z_Re_min`, `z_Re_max`, `z_Im_min`, `z_Im_max`: Define the calculation's complex plane boundaries.

## Performance

The program's parallel performance has been optimized by reducing unnecessary data transfer using point-to-point communication instead of collective calls (`MPI_Reduce`). This improves scalability for larger process counts.

## License

This project is open-source and free to use or modify for any purpose.
