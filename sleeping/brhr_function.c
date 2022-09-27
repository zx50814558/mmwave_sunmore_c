#include <stdlib.h>
#include <stdio.h> 
#include <stdint.h>
#include <math.h>
#include "brhr_function.h"

double *MLR_brhr(double *input, int delta, int len_)
{
	double *data_s = input;

	// Dynamic memory management
	double mean_ptr[len_];
	double m_ptr[len_];
	double b_ptr[len_];

	// Copy data from input data
	for (int i = 0; i < len_; i++)
	{
		mean_ptr[i] = input[i];
		m_ptr[i] = input[i];
		b_ptr[i] = input[i];
	}

	// Part1
	for (int t = 0; t < len_; t++)
	{
		if (!((t - delta) < 0 || (t + delta + 1) > len_))
		{
			int star = t - delta;
			int end = t + delta + 1;
			double mean_tmp = 0;
			for (int j = 0; j < end - star; j++)
				mean_tmp += input[star + j];
			mean_ptr[t] = (mean_tmp / (end - star));
			double mtmp = 0;
			for (int i = -delta; i < delta + 1; i++)
				mtmp += i * (input[t + i] - mean_ptr[t]);
			m_ptr[t] = (3 * mtmp) / (delta * (2 * delta + 1) * (delta + 1));
			b_ptr[t] = mean_ptr[t] - (t * m_ptr[t]);
		}
	}

	// Part2
	for (int t = 0; t < len_; t++)
	{
		if (!((t - delta) < 0 || (t + delta + 1) > len_))
		{
			double tmp_s = 0;
			for (int i = t - delta; i < t + delta; i++)
			{
				tmp_s += (m_ptr[i] * t) + b_ptr[i];
			}
			data_s[t] = tmp_s / (2 * delta + 1);
		}
	}

	return data_s;
}

static void lfilter(double *b, double *a, double *x, double *y, double *Z, int len_b, uint32_t len_x, int stride_X, int stride_Y)
{

	double *ptr_x = x, *ptr_y = y;
	double *ptr_Z;
	double *ptr_b = (double *)b;
	double *ptr_a = (double *)a;
	double *xn, *yn;
	const double a0 = *((double *)a);
	int n;
	uint32_t k;

	/* normalize the filter coefs only once. */
	for (n = 0; n < len_b; ++n)
	{
		ptr_b[n] /= a0;
		ptr_a[n] /= a0;
	}

	for (k = 0; k < len_x; k++)
	{
		ptr_b = (double *)b; /* Reset a and b pointers */
		ptr_a = (double *)a;
		xn = (double *)ptr_x;
		yn = (double *)ptr_y;
		if (len_b > 1)
		{
			ptr_Z = ((double *)Z);
			*yn = *ptr_Z + *ptr_b * *xn; /* Calculate first delay (output) */
			ptr_b++;
			ptr_a++;
			/* Fill in middle delays */
			for (n = 0; n < len_b - 2; n++)
			{
				*ptr_Z =
					ptr_Z[1] + *xn * (*ptr_b) - *yn * (*ptr_a);
				ptr_b++;
				ptr_a++;
				ptr_Z++;
			}
			/* Calculate last delay */
			*ptr_Z = *xn * (*ptr_b) - *yn * (*ptr_a);
		}
		else
		{
			*yn = *xn * (*ptr_b);
		}
		ptr_y += stride_Y; /* Move to next input/output point */
		ptr_x += stride_X;
	}
}

void _local_maxima_1d_brhr(double *x, int len_s, int *midpoints, int *m)
{

	// Dynamic memory management
	int left_edges[400] = {0};  // default 399
	int right_edges[400] = {0};  // default 399
	int i_ahead = 0;

	// Pointer
	int i = 1;			   // Pointer to current sample, first one can't be maxima
	int i_max = len_s - 1; // Last sample can't be maxima

	while (i < i_max)
	{

		// Test if previous sample is smaller
		if (x[i - 1] < x[i])
		{
			i_ahead = i + 1; // Index to look ahead of current sample

			// Find next sample that is unequal to x[i]
			while (i_ahead < i_max && x[i_ahead] == x[i])
				i_ahead += 1;

			// Maxima is found if next unequal sample is smaller than x[i]
			if (x[i_ahead] < x[i])
			{
				left_edges[*m] = i;
				right_edges[*m] = i_ahead - 1;
				midpoints[*m] = floor((left_edges[*m] + right_edges[*m]) / 2);
				// printf("midpoints[m]: %d\n", midpoints[m]);
				*m += 1;

				// Skip samples that can't be maximum
				i = i_ahead;
			}
		}
		i += 1;
	}
}

/* ----------- QuickSort (Start) ----------- */
// function to swap elements
void swap_brhr(int *a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

// function to find the partition position
int partition_brhr(int array[], int low, int high)
{

	// select the rightmost element as pivot
	int pivot = array[high];

	// pointer for greater element
	int i = (low - 1);

	// traverse each element of the array
	// compare them with the pivot
	for (int j = low; j < high; j++)
	{
		if (array[j] <= pivot)
		{

			// if element smaller than pivot is found
			// swap it with the greater element pointed by i
			i++;

			// swap element at i with element at j
			swap_brhr(&array[i], &array[j]);
		}
	}

	// swap the pivot element with the greater element at i
	swap_brhr(&array[i + 1], &array[high]);

	// return the partition point
	return (i + 1);
}

void quickSort_brhr(int array[], int low, int high)
{
	if (low < high)
	{

		// find the pivot element such that
		// elements smaller than pivot are on left of pivot
		// elements greater than pivot are on right of pivot
		int pi = partition_brhr(array, low, high);

		// recursive call on the left of pivot
		quickSort_brhr(array, low, pi - 1);

		// recursive call on the right of pivot
		quickSort_brhr(array, pi + 1, high);
	}
}

void brhr_function(double *sig, int input_len, int brhr, int *top, int *top_index)
{
    double forward, backward;
    double copy_sigs[input_len];
    copy_sigs[input_len - 1] = sig[input_len - 1];
    // --------------------- Phase_difference --------------------- 
    int len_phase_diff = input_len;
    for (int num = 1; num < len_phase_diff; num++)
        copy_sigs[num - 1] = sig[num] - sig[num - 1];
    
    // --------------------- RemoveImpulseNoise (RIN) ---------------------
    double RIN_par = 1.5;
    double removed_noise[input_len-1];
    removed_noise[0] = copy_sigs[0];
    for (int num = 1; num < len_phase_diff - 1; num++)
    {
        forward = copy_sigs[num] - copy_sigs[num - 1];
        backward = copy_sigs[num] - copy_sigs[num + 1];
        if ((forward > RIN_par && backward > RIN_par) || (forward < -RIN_par && backward < -RIN_par))
            removed_noise[num] = copy_sigs[num - 1] + (copy_sigs[num + 1] - copy_sigs[num - 1]) / 2;
        removed_noise[num] = copy_sigs[num];
    }

    // --------------------- iir_bandpass_filter_1 ---------------------
    //  if (order == 5) => BR 
    double y[input_len-1];
    if (brhr == 0) {
        double b[11] = {0.0003100856139325839, -0.0024503932074526414, 0.008191253474064738, -0.014462792713750464, 0.012602969766102465, 0.0, -0.012602969766102468, 0.01446279271375046, -0.008191253474064738, 0.0024503932074526414, -0.00031008561393258384};
        double a[11] = {1., -9.780812442849507, 43.08317719449593, -112.54898060868825, 193.10308878043222, -227.3654181511635, 186.0550507759317, -104.48315327333019, 38.53588426508279, -8.42919504975534, 0.8303585098573365};
        double delay[10] = {0}; // length of (a or b) - 1
        size_t len_ab = sizeof(a) / sizeof(double);
        double after_lifter[input_len-1];
        lfilter(b, a, removed_noise, y, delay, len_ab, input_len-1, 1, 1);
    }

    //  if (order == 9) => HR
    else {
        double b[19] = {0.0009309221423942934, -0.012651127859899214, 0.08140072903422219, -0.3279827818967048, 0.9206296690623369, -1.8881391347626006, 2.864155296830487, -3.1158825565644306, 2.0796297285426384, -1.6933375125414267e-15, -2.0796297285426384, 3.1158825565644315, -2.8641552968304875, 1.8881391347626006, -0.9206296690623373, 0.3279827818967048, -0.08140072903422219, 0.012651127859899216, -0.0009309221423942937};
        double a[19] = {1.0, -15.142612391789038, 109.52867216475022, -502.6626423702458, 1639.7914526180646, -4037.095860679349, 7772.407643496967, -11963.349399532222, 14922.878739349395, -15197.385490007382, 12665.576290591296, -8617.837327643654, 4751.9970631301, -2094.9246186588152, 722.2258948681153, -187.91289244221497, 34.75519411187606, -4.078772963228582, 0.22866640874033417};
        double delay[18] = {0};  // length of (a or b) - 1
        size_t len_ab = sizeof(a) / sizeof(double);
        double after_lifter[input_len-1];
        lfilter(b, a, removed_noise, y, delay, len_ab, input_len-1, 1, 1);
    }
    
    // --------------------- iir_bandpass_filter_1 ---------------------
    int iir_par = 2;
    int len_input = sizeof(y) / sizeof(double);
    double *data_s = MLR_brhr(y, iir_par, len_input);  // Output: data_s

    // --------------------- Feature_detection ---------------------
    // Signal length and half length
    int len_s_half = floor(len_input / 2);

    // Output peak
    int m_p = 0;  // Pointer to the end of valid area in allocated arrays
    int midpoints_peak[len_s_half];
    _local_maxima_1d_brhr(data_s, len_input, midpoints_peak, &m_p);  // Output: midpoints_peak

    // Dynamic memory management for negative signals
    double neg_x[len_input];
    for (int i = 0; i < len_input; i++)
        neg_x[i] = -data_s[i];

    // Output valley
    int m_v = 0;  // Pointer to the end of valid area in allocated arrays
    int midpoints_valley[len_s_half];
    _local_maxima_1d_brhr(neg_x, len_input, midpoints_valley, &m_v);  // Output: midpoints_valley

    // --------------------- Feature compress --------------------- 
    // Initialize (feature_compress)
    int total_feature[1300] = {0};
    int FC_par = 0;
    if (brhr == 0)
        FC_par = 22;
    else if (brhr == 1)
        FC_par = 5;

    // Given value
    for (int i = 0; i < m_p; i++)
        total_feature[i] = midpoints_peak[i];
    for (int i = 0; i < m_v; i++)
        total_feature[m_p+i] = midpoints_valley[i];

    // Perform quicksort on data
    quickSort_brhr(total_feature, 0, m_p + m_v - 1);
    
    // Compress parameter
    int ltera = 0;
    int sum_v = 0;
    int sum_p = 0;

    int start_feature, ltera_add, end_feature, location;
    int feature_compress_peak[m_p + m_v];
    int feature_compress_valley[m_p + m_v];
    // Algorithm
    while (ltera < m_p + m_v - 1)
    {
        // Record start at valley or peak (peak:0 valley:1)
        start_feature = 1;
        for (int i = 0; i < m_p; i++){
            if (midpoints_peak[i] == total_feature[ltera]){
                start_feature = 0;
                break;
            }
        }

        ltera_add = ltera;
        while (total_feature[ltera_add+1]-total_feature[ltera_add]<FC_par){
            // skip the feature which is too close
            ltera_add += 1;
            // break the loop if it is out of boundary
            if(ltera_add >= (m_p + m_v - 1))
                break;
        }

        // Record end at valley or peak (peak:0 valley:1)
        end_feature = 1;
        for (int i = 0; i < m_p; i++){
            if (midpoints_peak[i] == total_feature[ltera_add]){
                end_feature = 0;
                break;
            }
        }

        // If it is too close
        if (ltera != ltera_add){
            
            // situation1: began with valley end with valley
            if (start_feature == 1 && end_feature == 1){
                // using the lowest feature as represent
                location = ltera;
                for (int c = ltera; c < ltera_add; c++){
                    if (data_s[total_feature[c]] < data_s[total_feature[location]])
                        location = c;
                }
                feature_compress_valley[sum_v] = total_feature[location];
                sum_v += 1;
            }

            // situation2: began with valley end with peak
            else if (start_feature == 1 && end_feature == 0){
                feature_compress_valley[sum_v] = total_feature[ltera];
                feature_compress_peak[sum_p] = total_feature[ltera_add];
                sum_v += 1;
                sum_p += 1;
            }

            // situation3: began with peak end with valley
            else if (start_feature == 0 && end_feature == 1){
                feature_compress_valley[sum_v] = total_feature[ltera_add];
                feature_compress_peak[sum_p] = total_feature[ltera];
                sum_v += 1;
                sum_p += 1;
            }

            // situation4: began with peak end with peak
            else if (start_feature == 0 && end_feature == 0) {
                location = ltera;
                for (int c = ltera; c < ltera_add; c++){
                    if (data_s[total_feature[c]] > data_s[total_feature[location]])
                        location = c;
                }
                feature_compress_peak[sum_p] = total_feature[location];
                sum_p += 1;
            }
            ltera = ltera_add;
            
        }

        else {
            // It is normal featur point
            if (start_feature == 1){
                feature_compress_valley[sum_v] = total_feature[ltera];
                sum_v += 1;
            }
            else {
                if (sum_p < len_s_half){
                    feature_compress_peak[sum_p] = total_feature[ltera];
                    sum_p += 1;
                }
            }
        }
        ltera += 1;
    }
    
    // --------------------- Feature sort --------------------- 
    int compress_feature[sum_p + sum_v];
    // Given value
    for (int i = 0; i < sum_p; i++)
        compress_feature[i] = feature_compress_peak[i];
    for (int i = 0; i < sum_v; i++)
        compress_feature[sum_p+i] = feature_compress_valley[i];
    
    // Perform quicksort on data
    quickSort_brhr(compress_feature, 0, sum_p + sum_v - 1);
    
    // --------------------- Candidate search --------------------- 
    int CS_par;
    if (brhr == 0)
        CS_par = 17;
    else if (brhr == 1)
        CS_par = 4;

    double tmp_sum, window_sum, tmp_var, window_var;
    int NT_index, NB_index;
    double signal_pad[1300];
    int cadaidate_par = CS_par;
    // Doing the zero paddding
    for (int i = 0; i < len_input + (2 * cadaidate_par); i++) {
        if (i >= cadaidate_par && i < len_input + (2 * cadaidate_par) - cadaidate_par)
            signal_pad[i] = data_s[i-cadaidate_par];
        else
            signal_pad[i] = 1;
    }
    for (int i = 0; i < cadaidate_par; i++)
        signal_pad[i] = data_s[0];
    for (int i = len_input + (2 * cadaidate_par) - cadaidate_par; i < len_input + (2 * cadaidate_par) - 1; i++)
        signal_pad[i] = data_s[len_input-1];
    
    // Calaulate the mean and std using windows(for peaks)
    int NT_point[500];
    int NB_point[500];
    NB_index = 0;
    NT_index = 0;
    for (int i = 0; i < sum_v+sum_p; i++){
        // For the mean
        tmp_sum = 0;
        window_sum = 0;
        tmp_var = 0;
        window_var = 0;
        for (int j = compress_feature[i]; j < compress_feature[i]+2*cadaidate_par+1; j++)
            tmp_sum += signal_pad[j];
        window_sum = tmp_sum / (cadaidate_par*2+1);
        for (int j = compress_feature[i]; j < compress_feature[i]+2*cadaidate_par+1; j++)
            tmp_var += pow(signal_pad[j] - window_sum, 2) ;
        window_var = sqrt(tmp_var / (cadaidate_par*2+1));

        // determine if it is NT
        if (data_s[compress_feature[i]] > window_sum && window_var > 0.01){
            NT_point[NT_index] = compress_feature[i];
            NT_index += 1;
        }
        else if (data_s[compress_feature[i]] < window_sum && window_var > 0.01){
            NB_point[NB_index] = compress_feature[i];
            NB_index += 1;
        }
    }
    for (int i = 0; i < NT_index; i++) {
        top[i] = NT_point[i];
    }
    *top_index = NT_index;
}

void ada(double *sig, int input_len, int brhr, double *output)
{
    int top[500];
    int top_index = 0;
    int cur_index, pre_index;
    brhr_function(sig, input_len, brhr, top, &top_index);
    for (int i = 0; i < top_index - 1; i++) {
        cur_index = top[i + 1];
        pre_index = top[i];
        *output += fabs(sig[cur_index] - sig[pre_index]);
    }
}