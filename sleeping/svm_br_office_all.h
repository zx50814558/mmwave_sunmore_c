#include <stdlib.h>
#include <stdio.h> 
#include <math.h>

#define N_FEATURES_BR 3
#define N_CLASSES_BR 2
#define N_VECTORS_BR 166
#define N_ROWS_BR 2
#define N_COEFFICIENTS_BR 1
#define N_INTERCEPTS_BR 1
#define KERNEL_TYPE_BR 'r'
#define KERNEL_GAMMA_BR 0.01
#define KERNEL_COEF_BR 0.0
#define KERNEL_DEGREE_BR 3

int predict_br (double features_br[]);