
/******************************************************************
Description: Program to unit test the mandelbrot_mw.c using CUnit

Notes: compile with: gcc -o mandelbrot_test mandelbrot_mw.c -lm

******************************************************************/

#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "mandelbrot_mw.h"

/* Test configuration */
#define TEST_N_RE 120
#define TEST_N_IM 80

/* Test data */
int **test_nIter;
float *test_z_Re, *test_z_Im;

/* Helper function to initialize test data */
void setup_test_data()
{
    test_nIter = (int**)malloc(sizeof(int*) * (TEST_N_RE + 1));
    for (int i = 0; i < TEST_N_RE + 1; i++)
    {
        test_nIter[i] = (int*)calloc((TEST_N_IM + 1), sizeof(int));
    }
    test_z_Re = (float*)calloc((TEST_N_RE + 1), sizeof(float));
    test_z_Im = (float*)calloc((TEST_N_IM + 1), sizeof(float));

    calc_pnt_on_axs(test_z_Re, TEST_N_RE + 1, TEST_N_RE, z_Re_max, z_Re_min);
    calc_pnt_on_axs(test_z_Im, TEST_N_IM + 1, TEST_N_IM, z_Im_max, z_Im_min);
    // for (int i = 0; i < TEST_N_RE + 1; i++)
    // {
    //     test_z_Re[i] = (((float)i) / ((float)TEST_N_RE)) * (z_Re_max - z_Re_min) + z_Re_min;
    // }
    // for (int j = 0; j < TEST_N_IM + 1; j++)
    // {
    //     test_z_Im[j] = (((float)j) / ((float)TEST_N_IM)) * (z_Im_max - z_Im_min) + z_Im_min;
    // }
}

/* Test for calc_vals function */
void test_calc_vals()
{
    setup_test_data();
    /* Replace the global arrays in the original function with test data */
    for (int i = 0; i < TEST_N_RE + 1; i++)
    {
        calc_vals(i, TEST_N_IM + 1, test_nIter, test_z_Re, test_z_Im);
    }

    /* Verify that the iteration counts make sense */
    for (int i = 0; i < TEST_N_RE + 1; i++)
    {
        for (int j = 0; j < TEST_N_IM + 1; j++)
        {
            CU_ASSERT(test_nIter[i][j] >= 0);   // Iteration counts should not be negative
            CU_ASSERT(test_nIter[i][j] <= 100); // Iteration counts should not exceed the maximum
        }
    }
}

/* Test for write_to_file function */
void test_write_to_file()
{
    setup_test_data();

    // Call the write_to_file function
    write_to_file("TEST.dat", test_z_Re, test_z_Im, test_nIter, TEST_N_RE + 1, TEST_N_IM + 1);

    // Verify the output file contents
    FILE *file = fopen("TEST.dat", "r");
    CU_ASSERT_PTR_NOT_NULL(file); // Check that the file was created

    char line[256];
    int i = 0, j = 0;
    while (fgets(line, sizeof(line), file))
    {
        float re, im;
        int iter;

        // Parse the line
        int items = sscanf(line, "%f %f %d", &re, &im, &iter);
        CU_ASSERT_EQUAL(items, 3); // Check correct parsing

        // Verify the content matches the test data
        CU_ASSERT_DOUBLE_EQUAL(re, test_z_Re[i], 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(im, test_z_Im[j], 1e-6);
        CU_ASSERT_EQUAL(iter, test_nIter[i][j]);

        // Update indices
        j++;
        if (j > TEST_N_IM)
        {
            j = 0;
            i++;
        }
    }

    // Ensure we read the correct number of lines
    CU_ASSERT_EQUAL(i, TEST_N_RE + 1);

    fclose(file);

    // Clean up
    remove("TEST.dat");
}

/* Main function to set up and run tests */
int main(int argc, char *argv[])
{
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    CU_pSuite suite = CU_add_suite("Mandelbrot_Test_Suite", NULL, NULL);
    if (NULL == suite)
    {
        CU_cleanup_registry();
        return CU_get_error();
    }


    if ((NULL == CU_add_test(suite, "test of calc_vals", test_calc_vals)) ||
        (NULL == CU_add_test(suite, "test of write_to_file", test_write_to_file)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();

    return CU_get_error();
}
