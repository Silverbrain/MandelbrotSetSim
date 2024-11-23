/******************************************************************
Description: Program to calculate the Mandelbrot set using
			 a manager-worker pattern

Notes:
	compile with: mpicc -o mandelbrot_mw mandelbrot_mw.c -lm
	or simply run 'make' command in your terminal
******************************************************************/

#include <stdio.h>
#include <complex.h>
#include <stdbool.h>
#include <mpi.h>
#include "mandelbrot_mw.h"

//****************************/// WARNING ///**********************************
//	The following function is obsolete due to the new implementation method
//	and the `do_communication` should not be called. The code is being kept
//	only for backward compatibility.
//*****************************************************************************

/******************************************************************************
Communicate results so that the rank 1 process in the world group is ready
to write results out to file. A new communicator will be set up which includes
only the worker processes
******************************************************************************/
void do_communication(int myRank)
{

	MPI_Group worldGroup, workerGroup;
	MPI_Comm workerComm;
	int zeroArray = {0};

	int sendBuffer[(N_RE + 1) * (N_IM + 1)];
	int receiveBuffer[(N_RE + 1) * (N_IM + 1)];
	int index = 0;
	int i, j;

	// Get a group handle for the world group
	MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
	// Form a new group excluding world group rank 0
	MPI_Group_excl(worldGroup, 1, &zeroArray, &workerGroup);
	// Create a communicator for the new group
	MPI_Comm_create(MPI_COMM_WORLD, workerGroup, &workerComm);

	/* Pack nIter into a 1D buffer for sending*/
	for (i = 0; i < N_RE + 1; i++)
	{
		for (j = 0; j < N_IM + 1; j++)
		{
			sendBuffer[index] = nIter[i][j];
			index++;
		}
	}

	/* call MPI_reduce to collate all results on world group process 1
		The world group rank zero process does not make this call */
	if (myRank != 0)
	{
		MPI_Reduce(&sendBuffer, &receiveBuffer, (N_RE + 1) * (N_IM + 1), MPI_INT, MPI_SUM, 0, workerComm);
	}

	/* Unpack receive buffer into nIter */
	index = 0;
	for (i = 0; i < N_RE + 1; i++)
	{
		for (j = 0; j < N_IM + 1; j++)
		{
			nIter[i][j] = receiveBuffer[index];
			index++;
		}
	}

	/* Free the group and communicator */
	if (myRank != 0)
	{
		MPI_Comm_free(&workerComm);
	}
	MPI_Group_free(&workerGroup);
}

/******************************************************************************************/

int main(int argc, char *argv[])
{

	/* MPI related variables */
	int myRank;			 /* Rank of this MPI process */
	int nProcs;			 /* Total number of MPI processes*/
	int nextProc;		 /* Next process to send work to */
	MPI_Status status;	 /* Status from MPI calls */
	int endFlag = -9999; /* Flag to indicate completion*/

	/* Timing variables */
	double start_time, end_time;

	/* Loop indices */
	int i, j;

	MPI_Init(&argc, &argv);

	/* Record start time */
	start_time = MPI_Wtime();

	/* Get job size and rank information */
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	if (myRank == 0 && verbose)
	{
		printf("Calculating Mandelbrot set with %d processes\n", nProcs);
	}

	/* Initialise nIter with zero to collate results into it */
	init_arr();

	calc_pnt_on_axs(z_Re, N_RE + 1, (float)N_RE, z_Re_max, z_Re_min); // Points on real axis
	calc_pnt_on_axs(z_Im, N_IM + 1, (float)N_IM, z_Im_max, z_Im_min); // Points on imaginary axis
	//calc_pnt_on_IM_axs(); // Points on imaginary axis

	// Message buffer for Manager/Worker communications
	int res_buff[N_IM + 3]; // calculated column

	// Constants pointing to the values encapsulated in the `res_buffer`
	const int col_idx = N_IM + 1;	// index of the column calculated by worker
	const int proc_rank = N_IM + 2; // Rank of the worker that calculated the column

	// Initialise the message buffer
	res_buff[col_idx] = -1; // Set the initial value of column index to an invalid value
	res_buff[proc_rank] = myRank;

	for (i = 0; i < N_IM + 1; i++)
	{
		res_buff[i] = 0;
	}

	// Manager process
	if (myRank == 0)
	{
		// Hand out work to worker processes
		for (i = 0; i < N_RE + 1; i++)
		{
			// Receive request for work
			MPI_Recv(&res_buff, N_IM + 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Extracting the rank of nextProc and index of computed column from message buffer*/
			unpack_msg(&nextProc, res_buff, proc_rank, col_idx);

			// Send i value to requesting process
			MPI_Send(&i, 1, MPI_INT, nextProc, 100, MPI_COMM_WORLD);
		}
		// Tell all the worker processes to finish (once for each worker process)
		for (i = 0; i < nProcs - 1; i++)
		{
			// Receive request for work
			MPI_Recv(&res_buff, N_IM + 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/* Extracting the rank of nextProc and index of computed column from message buffer*/
			unpack_msg(&nextProc, res_buff, proc_rank, col_idx);
			
			// Send endFlag to finish
			MPI_Send(&endFlag, 1, MPI_INT, nextProc, 100, MPI_COMM_WORLD);
		}
	} // if manager process

	// Worker Processes
	else
	{
		while (true)
		{
			// Send the computed column and request for work
			MPI_Send(&res_buff, N_IM + 3, MPI_INT, 0, 100 + myRank, MPI_COMM_WORLD);
			// Receive i value to work on
			MPI_Recv(&i, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);

			if (i == endFlag)
			{
				break;
			}
			else
			{
				// Set the index of the column which is being computed to send to the Manager process
				res_buff[col_idx] = i;

				calc_vals(i, N_IM + 1, nIter, z_Re, z_Im);

				// pack the calculated column into the message buffer
				for (j = 0; j < N_IM + 1; j++)
				{
					res_buff[j] = nIter[i][j];
				}
			}

		} // while(true)
	} // else worker process

	/*
		Since the communication pattern has changed,
		there is no need to call `do_communication` function,
		anymore therefore call to the
		`do_communication` function is commented out.
	*/
	// do_communication(myRank); // Communicate results so rank 0 has all the values

	/* Write out results */
	if (doIO && myRank == 0)
	{
		if (verbose)
		{
			printf("Writing out results from process %d \n", myRank);
		}
		write_to_file("mandelbrot_sim.dat", z_Re, z_Im, nIter, N_RE + 1, N_IM + 1);
	}

	/* Record end time */
	MPI_Barrier(MPI_COMM_WORLD);
	end_time = MPI_Wtime();

	/* Record end time. The barrier synchronises the process so they all measure the same time */
	if (myRank == 0)
	{
		printf("STATS (num procs, elapsed time): %d %f\n", nProcs, end_time - start_time);
	}

	// Cleanup memory
	clean_arr();

	MPI_Finalize();

	return 0;
}
