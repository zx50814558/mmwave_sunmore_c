#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* 
Used to calculate respiratory rate and heart rhythm
sig = input
input_len = input data size
brhr = Switching breathing(0) or heartbeat(1)
top = output top feature
top_index = output number of top feature in top
return = None
*/
void brhr_function(double *sig, int input_len, int brhr, int *top, int *top_index);

/*
Used to calculate ada
sig = input
input_len = input data size
brhr = Switching breathing(0) or heartbeat(1)
output = ada output
return = None
*/
void ada(double *sig, int input_len, int brhr, double *output);