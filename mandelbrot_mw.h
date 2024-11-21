#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

/* Number of intervals on real and imaginary axes*/
#define N_RE 12000
#define N_IM 8000

/* Number of iterations at each z value */
int **nIter;

/* Points on real and imaginary axes*/
float *z_Re, *z_Im;

/* Domain size */
const float z_Re_min = -2.0; /* Minimum real value*/
const float z_Re_max = 1.0;  /* Maximum real value */
const float z_Im_min = -1.0; /* Minimum imaginary value */
const float z_Im_max = 1.0;  /* Maximum imaginary value */

/* Set to true to write out results*/
const bool doIO = true;
const bool verbose = true;

/******************************************************************************************
    Calculate number of iterations for a given `i` index (real axis) and loop over all j values
    (imaginary axis specified by `col_size`). The `nIter` array is populated with the results.
******************************************************************************************/
void calc_vals(int i, int col_size, int **_nIter, float *_z_Re, float *_z_Im)
{
    /* Maximum number of iterations*/
    const int maxIter = 100;

    /* Value of Z at current iteration*/
    float complex z;
    /*Value of z at iteration zero*/
    float complex z0;

    int j, k;

    /* Loop over imaginary axis */
    for (j = 0; j < col_size; j++)
    {
        z0 = _z_Re[i] + _z_Im[j] * I;
        z = z0;

        /* Iterate up to a maximum number or bail out if mod(z) > 2 */
        k = 0;
        while (k < maxIter)
        {
            _nIter[i][j] = k;
            if (cabs(z) > 2.0)
                break;
            z = z * z + z0;
            k++;
        }
    }
}

/******************************************************************************************
Write the results of the simulation to the specified `filename`.
******************************************************************************************/
void write_to_file(char filename[], float _z_Re[], float _z_Im[], int *_nIter[], int n, int m)
{
    int i, j;
    FILE *outfile;

    outfile = fopen(filename, "w");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            fprintf(outfile, "%f %f %d \n", _z_Re[i], _z_Im[j], _nIter[i][j]);
        }
    }
    fclose(outfile);
}

/******************************************************************************************
Initialise the input array with zero to avoid arbiterary values
******************************************************************************************/
void init_arr()
{
    nIter = (int **)malloc(sizeof(int *) * (N_RE + 1));
    for (int i = 0; i < N_RE + 1; i++)
    {
        nIter[i] = (int *)calloc((N_IM + 1), sizeof(int));
    }
    z_Re = (float *)calloc((N_RE + 1), sizeof(float));
    z_Im = (float *)calloc((N_IM + 1), sizeof(float));
}

/******************************************************************************************
Free up the allocated memory used in the proram.
******************************************************************************************/
void clean_arr()
{
    // Free the memory
    for (int i = 0; i < N_RE + 1; i++)
    {
        free(nIter[i]); // Free each row
    }
    free(nIter);

    free(z_Re);
    free(z_Im);
}

/******************************************************************************************
Calculate points on axis
******************************************************************************************/
void calc_pnt_on_axs(float *z, int size, float N, float z_max, float z_min)
{
    for (int i = 0; i < size; i++)
    {
        z[i] = (((float)i) / N) * (z_max - z_min) + z_min;
    }
}

/******************************************************************************************
Unpack recieved message to extract calculated column and next process
******************************************************************************************/
void unpack_msg(int *nextProc, int *res_buff, int proc_rank, int col_idx)
{
    /* Extracting the rank of nextProc and index of computed column from message buffer*/
    *nextProc = res_buff[proc_rank];
    int calc_col = res_buff[col_idx];

    // Check if the column index is valid
    if (calc_col >= 0)
    {
        // Unpack the message buffer and store the received column in the nIter array
        for (int j = 0; j < N_IM + 1; j++)
        {
            nIter[calc_col][j] = res_buff[j];
        }
    }
}
