#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N_FEATURES_HR 3
#define N_CLASSES_HR 2
#define N_VECTORS_HR 173
#define N_ROWS_HR 2
#define N_COEFFICIENTS_HR 1
#define N_INTERCEPTS_HR 1
#define KERNEL_TYPE_HR 'r'
#define KERNEL_GAMMA_HR 0.1
#define KERNEL_COEF_HR 0.0
#define KERNEL_DEGREE_HR 3

int predict_hr (double features_hr[]);