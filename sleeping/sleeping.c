//  C library headers
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

// Linux headers
#include <fcntl.h>	 // Contains file controls like O_RDWR
#include <errno.h>	 // Error integer and strerror() function
#include <termios.h> // Contains POSIX terminal control definitions
#include <unistd.h>	 // write(), read(), close()
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <complex.h>

// Additional function
#include "pocketfft.h"
#include "polyfit.h"
#include "brhr_function.h"

// Sklearn model
#include "svm_br_office_all.h"
#include "svm_hr_office_all.h"
#include "sleep_feature_min_rf.h"

// Helper functions to get the min and max of two numbers
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))

struct timeval start, stop;

// Data containers
double unwrapPhasePeak_mm[800];
double heartRateEst_FFT_mean[800];
double heartRateEst_xCorr_mean[800];
double breathingEst_FFT_mean[800];
double breathingEst_xCorr_mean[800];
double breath_ti[800];
double heart_ti[800];
double current_window_bmi[1200];  // raw_sig[-60*20:]
double LF_HF_LFHF_windows[6100];  // raw_sig[-5*60*20:]

// Sklearn to c
int svm_result;
double svm_input[3];


double *MLR(double *input, int delta, int len_)
{
	double *data_s = input;

	// Dynamic memory management
	double mean_ptr[799] = {0};
	double m_ptr[799] = {0};
	double b_ptr[799] = {0};

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

int *_local_maxima_1d(double *x, int len_s, int *midpoints, int *m)
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
	return midpoints;
}

/* ----------- QuickSort (Start) ----------- */
// function to swap elements
void swap(int *a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

// function to find the partition position
int partition(int array[], int low, int high)
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
			swap(&array[i], &array[j]);
		}
	}

	// swap the pivot element with the greater element at i
	swap(&array[i + 1], &array[high]);

	// return the partition point
	return (i + 1);
}

void quickSort(int array[], int low, int high)
{
	if (low < high)
	{

		// find the pivot element such that
		// elements smaller than pivot are on left of pivot
		// elements greater than pivot are on right of pivot
		int pi = partition(array, low, high);

		// recursive call on the left of pivot
		quickSort(array, low, pi - 1);

		// recursive call on the right of pivot
		quickSort(array, pi + 1, high);
	}
}
/* ----------- QuickSort (End) ----------- */

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

typedef union
{
	float f;
	struct
	{
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} raw;
} myfloat;

unsigned int convertToInt(int *arr, int low, int high)
{
	unsigned f = 0, i;
	for (i = high; i >= low; i--)
	{
		f = f + arr[i] * pow(2, high - i);
	}
	return f;
}

unsigned int ieee754_convert(int *ieee)
{
	myfloat var;
	unsigned f = convertToInt(ieee, 9, 31);
	var.raw.mantissa = f;
	f = convertToInt(ieee, 1, 8);
	var.raw.exponent = f;
	var.raw.sign = ieee[0];
	return var.f;
}

void array_shift()
{
	for (int num = 0; num < 799; num++)
	{
		unwrapPhasePeak_mm[num] = unwrapPhasePeak_mm[num + 1];
		heartRateEst_FFT_mean[num] = heartRateEst_FFT_mean[num + 1];
		heartRateEst_xCorr_mean[num] = heartRateEst_xCorr_mean[num + 1];
		breathingEst_FFT_mean[num] = breathingEst_FFT_mean[num + 1];
		breathingEst_xCorr_mean[num] = breathingEst_xCorr_mean[num + 1];
		breath_ti[num] = breath_ti[num + 1];
		heart_ti[num] = heart_ti[num + 1];
	}
}

// 當呼吸律或心律為異常值, 以前一秒的輸出值取代。1 = BR 0 = HR
double substitute(double pre, double input, int result_type){
	if (result_type == 0){
		if (input < 40 || input > 110)
			return pre;
		else
			return input;
	}
	else {
		if (input < 10 || input > 25)
			return pre;
		else
			return input;
	}
}

// Sleeping features function
void mov_dens_fn(double *raw_sig, double *percent){
    int x_index;
    double tmp_x, tmp_top, x_mean, result;
    double count = 0;
    double x[4] = {0};
    double top[4] = {0};
    
    for (int num = 0; num < 120; num++){
        x_mean = 0;
        x_index = 0;
        result = 0;
        for (int index = num*4; index < (num+1)*4; index++){
            tmp_x = round(raw_sig[index] * 100000000) / 100000000;
            x[x_index] = tmp_x;
            x_mean += tmp_x;
            x_index++;
        }
        x_mean = x_mean / 4;
        for (int i = 0; i < 4; i++){
            tmp_top = pow(x[i] - x_mean, 2);
            top[i] = tmp_top;
            result += tmp_top;
        }
        result = result / 3;
        if (result > 0.045)
            count++;
    }
    *percent = (count / 120) * 100;
}

// tfRSA_fn(intput, output)
//> intput = array
//> output = features
void tfRSA_fn(double *fRSA_sig, double *tfRSA){
    double sum = 0.0, mean, SD = 0.0;
    for (int i = 0; i < 10; ++i) {
        sum += fRSA_sig[i];
    }
    mean = sum / 10;
    for (int i = 0; i < 10; ++i) {
        SD += pow(fRSA_sig[i] - mean, 2);
    }
    *tfRSA = sqrt(SD / 10);
}

void convolve(double *h, double *x, int lenH, int lenX, double *output)
{
    int nconv = lenH+lenX-1;
    int h_start, x_start, x_end;
    double y[100] = {0};
    for (int i = 0; i < nconv; i++){
        x_start = MAX(0, i - lenH + 1);
        x_end   = MIN(i + 1, lenX);
        h_start = MIN(i, lenH - 1);
        for (int j = x_start; j < x_end; j++)
            y[i] += h[h_start--] * x[j];
    }
    int start_index = (lenX-1) / 2;
    for (int i = 0; i < lenX; i++)
        output[i] = y[i + start_index];
}

// Only supple window_length = 31 and polyorder == (2 or 3)
void savgol_filter(double *x, int window_length, int polyorder, double *input_mean, double *y) {
    double conv_y[100] = {0};
    double coeffs_2[] = {-0.041055718475071855, -0.026392961876831954, -0.012741429871574048, -0.00010112245929822615, 0.011527960359995528, 0.02214581858630722, 0.031752452219636844, 0.04034786125998441, 0.0479320457073499, 0.05450500556173333, 0.060066740823134686, 0.06461725149155399, 0.06815653756699122, 0.07068459904944638, 0.07220143593891949, 0.0727070482354105, 0.07220143593891949, 0.07068459904944638, 0.06815653756699123, 0.06461725149155399, 0.06006674082313469, 0.05450500556173333, 0.04793204570734991, 0.040347861259984415, 0.03175245221963685, 0.022145818586307226, 0.01152796035999553, -0.00010112245929822268, -0.012741429871574055, -0.02639296187683194, -0.04105571847507189};
    double coeffs_3[] = {-0.04105571847507182, -0.02639296187683192, -0.012741429871574027, -0.00010112245929820307, 0.011527960359995546, 0.022145818586307233, 0.03175245221963686, 0.040347861259984415, 0.0479320457073499, 0.054505005561733336, 0.060066740823134686, 0.06461725149155399, 0.06815653756699121, 0.07068459904944636, 0.07220143593891948, 0.07270704823541049, 0.07220143593891948, 0.07068459904944638, 0.06815653756699121, 0.06461725149155398, 0.06006674082313469, 0.05450500556173333, 0.04793204570734991, 0.040347861259984415, 0.03175245221963686, 0.022145818586307237, 0.01152796035999555, -0.0001011224592981996, -0.012741429871574027, -0.026392961876831905, -0.041055718475071855};
    
    // convolve
    if (polyorder == 2){
        convolve(coeffs_2, x, window_length, window_length, conv_y);
    }
               
    else if (polyorder == 3){
        convolve(coeffs_3, x, window_length, window_length, conv_y);
    }

    // savgol_filter
    int rVal;
    double mean_tmp = 0;
    double poly_coeffs[5] = {0};
    double poly_x[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };

    rVal = polyfit( window_length, poly_x, x, polyorder+1, poly_coeffs);

    // np.polyval
    double input_first[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    double input_last[] = {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

    size_t input_size = sizeof(input_first) / sizeof(double);

    for (int i = 0; i < input_size; i++) {
        if (polyorder == 2)
            *(y+i) = pow(input_first[i], 2)*poly_coeffs[0] + pow(input_first[i], 1)*poly_coeffs[1] + pow(input_first[i], 0)*poly_coeffs[2];
        else if(polyorder == 3)
            *(y+i) = pow(input_first[i], 3)*poly_coeffs[0] + pow(input_first[i], 2)*poly_coeffs[1] + pow(input_first[i], 1)*poly_coeffs[2] + pow(input_first[i], 0)*poly_coeffs[3];
        *input_mean  += y[i];
    }
    *(y+15) = conv_y[15];
    *input_mean  += conv_y[15];
    for (int i = 0; i < input_size; i++) {
        if (polyorder == 2)
            *(y+i+input_size+1) = pow(input_last[i], 2)*poly_coeffs[0] + pow(input_last[i], 1)*poly_coeffs[1] + pow(input_last[i], 0)*poly_coeffs[2];
        else if (polyorder == 3)
            *(y+i+input_size+1) = pow(input_last[i], 3)*poly_coeffs[0] + pow(input_last[i], 2)*poly_coeffs[1] + pow(input_last[i], 1)*poly_coeffs[2] + pow(input_last[i], 0)*poly_coeffs[3];
        *input_mean  += y[i+input_size+1];
    }
    *input_mean = *input_mean  / 31;
}

void var_RPM(double *sig, double *output) {
    double rk_mean = 0;
    for (int i = 0; i < 600; i++)
        rk_mean += sig[i];
    rk_mean = rk_mean / 600;

    double tmp_mean;
    double tmp_output = 0;
    double brhr_mins[10] = {0};
    for (int i = 0; i < 10; i++) {
        tmp_mean = 0;
        for (int j = 0; j < 60; j++)
            tmp_mean += sig[(60*i)+j];
        tmp_output += pow((tmp_mean / 60) - 1 * rk_mean, 2);
    }
    *output = tmp_output / 9;
}

void bmi(double *sig, double *output) {
    double ak_array[6] = {0};
    
    double ak_min = 0;
    for (int i = 0; i < 10; i++)
        ak_min += sig[i];
    ak_min = ak_min / 10;

    double ak;
    for (int i = 0; i < 6; i++) {
        ak = 0;
        for (int j = 0; j < 10; j++) {
            ak += sig[(i*10)+j];
        }
        ak = (ak / 10);
        ak_array[i] = ak;
        if (ak < ak_min)
            ak_min = ak;
    }

    double result = 0;
    for (int i = 0; i < 6; i++)
        result += ak_array[i] - ak_min;
    *output = result;
}

void deep_parameter(double bi, double h, double *output) {
    *output = bi / (h + bi);
}

void time_fn (int hours, int minutes, int seconds, int *time_featres) {
    *time_featres = (hours + 24 - 20) * 3600 + minutes * 60 + seconds;
}

void rem_parameter(double *sig, double *output)
{
    double former;
    double latter;
    double total_avg = 0;
    for (int i = 0; i < 5; i++) {
        former = 0;
        latter = 0;
        for (int j = i*30; j < (i+1)*30; j++)
            former += sig[j];
        for (int j = (i+1)*30; j < (i+2)*30; j++)
            latter += sig[j];
        total_avg += fabs(former/30 - latter/30);
    }
    *output = total_avg / 5;
}

void mean_fn_double (double *sig, int sig_len, double *output_d)
{
    *output_d = 0;
    for (int i = 0; i < sig_len; i++)
        *output_d += sig[i];
    *output_d = *output_d / sig_len;
}

void mean_fn_int (int *sig, int sig_len, double *output_i)
{
    *output_i = 0;
    for (int i = 0; i < sig_len; i++)
        *output_i += sig[i];
    *output_i = *output_i / sig_len;
}

int main(void)
{
    char filename[100] = {0};
    char input_name[100];
    char *input_n;
    char *root_dir = "dataset/"; // 改檔名
    printf("Input file name = ");
    input_n = fgets(input_name, 100, stdin);
    input_n[strcspn(input_n, "\r\n")] = 0;
    strcat(filename, root_dir);
    strcat(filename, input_name);
    strcat(filename, ".csv");
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("error");
        return -1;
    }
	fprintf(fp, "heart, breath, bmi, deep_p, ada_br, ada_hr, var_RPM, var_HPM, rem_parameter, mov_dens, LF, HF, LFHF, sHF, sLFHF, tfRSA, tmHR, sfRSA, smHR, sdfRSA, sdmHR, stfRSA, stmHR, time, datetime, sleep\n");
	fclose(fp);
	int len_s_half;
	int *feature_peak, *feature_valley;
	int m_p, m_v; // Number of peak & valley
	double neg_x[799] = {0};

	// Initialize (feature_compress)
	int start_feature, end_feature, ltera_add, time_thr;
	int location, ltera, sum_v, sum_p;

    int serial_port = open("/dev/ttyTHS1", O_RDWR);
	struct termios tty;

    if (tcgetattr(serial_port, &tty) != 0)
	{
		printf("Error %i from tcgetattr: %s\n", errno, strerror(errno));
		return 1;
	}
	tty.c_cflag &= ~PARENB;
	tty.c_cflag &= ~CSTOPB;
	tty.c_cflag &= ~CSIZE;
	tty.c_cflag |= CS8;
	tty.c_cflag &= ~CRTSCTS;
	tty.c_cflag |= CREAD | CLOCAL;
	tty.c_lflag &= ~ICANON;
	tty.c_lflag &= ~ECHO;
	tty.c_lflag &= ~ECHOE;
	tty.c_lflag &= ~ECHONL;
	tty.c_lflag &= ~ISIG;
	tty.c_iflag &= ~(IXON | IXOFF | IXANY);
	tty.c_iflag &= ~(IGNBRK | BRKINT | PARMRK | ISTRIP | INLCR | IGNCR | ICRNL);
	tty.c_oflag &= ~OPOST;
	tty.c_oflag &= ~ONLCR;
	tty.c_cc[VTIME] = 10;
	tty.c_cc[VMIN] = 0;

	cfsetispeed(&tty, B921600);
	cfsetospeed(&tty, B921600);

	if (tcsetattr(serial_port, TCSANOW, &tty) != 0)
	{
		printf("Error %i from tcsetattr: %s\n", errno, strerror(errno));
		return 1;
	}

    char read_buf[1024];
	int size;
	int data_idx = 0;
	unsigned long int header_reader_output[10];
	unsigned long int byte[4];
	unsigned long int tlv_header_reader_output[2];
	unsigned int vsos_byte[128];
	unsigned int tlv_header[4];
	short int rangeProfile[2];
	float vsos_array[34];			   //主要輸出1
	short int rangeProfile_array[126]; //主要輸出2
	short int int16_number;
	int i = 0;
	int j = 0;
	int k = 0;
	int int16_temp_number = 0;
	unsigned int ieee[32];
	unsigned int int16_temp[16];
	int magicWord[8] = {2, 1, 4, 3, 6, 5, 8, 7};

	//讀值之後的變數
	time_t start_time;
	start_time = time(NULL);
	int array_index = 0;
	int array_index_bmi = 0;
	float phase_diff[799];
	double removed_noise[799];
	float forward, backward;
	float hr_mean_FFT, hr_mean_xCorr, br_mean_FFT, br_mean_xCorr, breath_mean_ti, heart_mean_ti;
	float thr;

	// 動態記憶體
	int total_feature[800] = {0};
	int feature_compress_valley[800] = {0};
	int feature_compress_peak[800] = {0};
	int compress_feature[800] = {0};
	int NT_point[800] = {0};
	int NB_point[800] = {0};
	double signal_pad[1600] = {0};

	// 最終生理資訊
	double final_hr, final_br, final_hr_sub, final_br_sub, br_rpm, hr_rpm;
	double tmp_br = 0;
	double tmp_hr = 0;

	// Time data
	int hours, minutes, seconds, day, month, year;
	time_t now;
	int start_time_sec, start_day, start_month, start_year;
	int next_YMD = 0;

	// Start record feature
	int counter = 0;
	int begin = 0;
	
	// Features (var_RPM)
	int var_index = 0;
	double var_RPM_br_KNN[10*60] = {0};
	double var_RPM_hr_KNN[10*60] = {0};

	// Features (mov_dens)
	int mov_dens_index = 0;
	double mov_dens = 0;
	double mov_dens_ar[60] = {0};

	// Features (tfRSA)
	int tfRSA_index = 0;
	int open_stfRSA = 0;
	double tfRSA = 0;
	double tfRSA_ar[31] = {0};
	double tfRSA_contain10_br[10] = {0};  // Contain 10 breath rates

	// Features (sfRSA)
	int sfRSA_index = 0;
	double sfRSA[31] = {0};
	double sfRSA_mean = 0;
	double sfRSA_contain31_br[31] = {0};  // Contain 31 breath rates
	double sfRSA_ar[60] = {0};

	// Features (stfRSA)
	int stfRSA_index = 0;
	double stfRSA[31] = {0};
	double stfRSA_mean = 0;
	double stfRSA_ar[60] = {0};

	// Features (sdfRSA)
	int sdfRSA_index = 0;
	double sdfRSA[31] = {0};
	double sdfRSA_mean = 0;
	double sdfRSA_tmp[31] = {0};
	double sdfRSA_ar[60] = {0};

	// Features (tmHR)
	int tmHR_index = 0;
	int open_stmHR = 0;
	double tmHR = 0;
	double tmHR_ar[31] = {0};
	double tmHR_contain10_hr[10] = {0};  // Contain 10 heart rates

	// Features (smHR)
	int smHR_index = 0;
	double smHR_contain31_hr[31] = {0};  // Contain 31 heart rates
	double smHR_mean = 0;
	double smHR[31] = {0};
	double smHR_ar[60] = {0};

	// Features (sdmHR)
	int sdmHR_index = 0;
	double sdmHR_tmp[31] = {0};
	double sdmHR_mean = 0;
	double sdmHR[31] = {0};
	double sdmHR_ar[60] = {0};

	// Features (stmHR)
	int stmHR_index = 0;
	double stmHR[31] = {0};
	double stmHR_mean = 0;
	double stmHR_ar[60] = {0};

	// Features (LF_HF_LFHF)
	int LF_HF_LFHF_index = 0;
	int sF_index = 0;
	double input_energe[6100] = {0};
	double out_fft[6100] = {0};
	double emerge_sum_LF = 0;
	double emerge_sum_HF = 0;
	double LFHF_eng = 0;
	double LF_ar[60] = {0};
	double HF_ar[60] = {0};
	double LFHF_ar[60] = {0};
	double HF_arr[31] = {0};
	double LFHF_arr[31] = {0};
	int open_HF = 0;
	double sHF[31] = {0};
	double sHF_mean = 0;
	double sLFHF[31] = {0};
	double sLFHF_mean = 0;
	double sHF_ar[60] = {0};
	double sLFHF_ar[60] = {0};
	int sHF_index = 0;

	// 睡眠階段 (Paper2)
	// Features (Variance of RPM)
	int var_RPM_index = 0;
	double var_RPM_br = 0;
	double var_RPM_hr = 0;
	double var_RPM_ar[60] = {0};
	double var_HPM_ar[60] = {0};
	int looper = 0;

	// Body Movement Index (BMI)
    // Features (Deep Parameter)
	double bmi_current = 0;
	double hk = 0;
	double dk = 0;
	double bmi_ar[60] = {0};
	double deep_p_ar[60] = {0};
	int bmi_ar_index = 0;

	// Features (ADA)
	int ada_index = 0;
	double ada_br_ar[60] = {0};
	double ada_hr_ar[60] = {0};

	// Features (REM_par)
	double rem_par = 0;
	double rem_contain[5*60] = {0};
	double rem_parameter_ar[60] = {0};
	int rem_par_index = 0;

	// Features (time_fn)
	int hours_tf, minutes_tf, seconds_tf, time_featres;
	time_t now_tf;
	int time_ar[60] = {0};
	int time_fn_index = 0;

	// Features (br_hr)
	double breath_ar[60] = {0};
	double heart_ar[60] = {0};
	int br_hr_index = 0;

	// time par
	int start_hour, start_min, end_hour, end_min;
	int counter_mean = 0;
	int next_HM = 0;
	
	// all_results
	double all_results[24] = {0};
	int predict_result;

	time(&now);
	struct tm *local = localtime(&now);
	hours = local->tm_hour;         // 獲取自午夜以來的小時數 (0-23)
    minutes = local->tm_min;        // 獲取小時後經過的分鐘數 (0-59)
    seconds = local->tm_sec;        // 獲取一分鐘後經過的秒數 (0-59)
	start_time_sec = seconds + minutes * 60 + hours * 3600;
    start_day = local->tm_mday;            // 獲取月份中的日期(1 到 31)
    start_month = local->tm_mon + 1;      // 獲取一年中的月份(0 到 11)
    start_year = local->tm_year + 1900;   // 獲取自 1900 年以來的年份
	

	while (1)
	{
		while ((size = read(serial_port, &read_buf, sizeof(read_buf) - 1)) > 0)
		{
			// printf("----------start-----------\n");
			gettimeofday(&start, NULL);
			int state = 0;
			if (state == 0)
			{
				int same = 0;
				for (int ix = 0; ix < 8; ++ix)
				{
					// printf("%d\n", read_buf[ix]);
					if (read_buf[ix] == magicWord[ix])
					{
						same++;
						if (same == 8)
						{
							state = 1;
						}
					}
				}
			}
			if (state == 1)
			{
				for (int ix = 0; ix < 40; ++ix)
				{
					byte[i] = read_buf[ix + 8];
					i++;
					if (i == 4)
					{
						header_reader_output[((ix + 1) / 4) - 1] = ((byte[0] << 0) + (byte[1] << 8) + (byte[2] << 16) + (byte[3] << 24));
						i = 0;
					}
				}
				for (int ix = 48; ix < 176; ++ix)
				{
					vsos_byte[j] = read_buf[ix];
					j++;
					if (j == 4)
					{
						for (int ix2 = 0; ix2 < 2; ix2++)
						{
							vsos_array[ix2] = (float)((vsos_byte[ix2 * 2] << 0) + (vsos_byte[(ix2 * 2) + 1] << 8));
						}
					}
					if (j == 8)
					{
						for (int ix1 = 7; ix1 > 3; --ix1)
						{
							for (int ix = 0; ix < 8; ++ix)
							{
								ieee[k] = (vsos_byte[ix1] << ix & 0x80) >> 7;
								k++;
							}
						}
						k = 0;
						float ieee775544 = ieee754_convert(ieee);
						vsos_array[2] = ieee775544;
					}
					if (j == 18)
					{
						for (int ix2 = 4; ix2 < 8; ix2++)
						{
							vsos_array[ix2 - 1] = (float)((vsos_byte[ix2 * 2] << 0) + (vsos_byte[(ix2 * 2) + 1] << 8));
						}
					}
					if (j == 128)
					{
						for (int ix2 = 0; ix2 < 27; ix2++)
						{
							for (int ix1 = 19 + (ix2 * 4); ix1 > 15 + (ix2 * 4); --ix1)
							{
								for (int ix = 0; ix < 8; ++ix)
								{
									ieee[k] = (vsos_byte[ix1] << ix & 0x80) >> 7;
									k++;
								}
							}
							k = 0;
							myfloat var;
							unsigned f = convertToInt(ieee, 9, 31);
							var.raw.mantissa = f;
							f = convertToInt(ieee, 1, 8);
							var.raw.exponent = f;
							var.raw.sign = ieee[0];
							vsos_array[ix2 + 7] = var.f; // 7-33
						}
					}
				}
				j = 0;
				for (int ix = 176; ix < 184; ++ix)
				{
					tlv_header[i] = read_buf[ix];
					i++;
					if (i == 4)
					{
						i = 0;
					}
				}

				for (int ix = 184; ix < 436; ++ix)
				{
					rangeProfile[i] = read_buf[ix];
					i++;
					if (i == 2)
					{
						int16_number = (rangeProfile[0] << 0) + (rangeProfile[1] << 8);
						rangeProfile_array[(ix - 185) / 2] = int16_number;
						i = 0;
					}
				}
			}
			break;
		}

		// Reading data
		time_t end_time;
		end_time = time(NULL);
		if (array_index < 800) {
			unwrapPhasePeak_mm[array_index] = (double)vsos_array[7];
			heartRateEst_FFT_mean[array_index] = (double)vsos_array[10];
			heartRateEst_xCorr_mean[array_index] = (double)vsos_array[12];
			breathingEst_FFT_mean[array_index] = (double)vsos_array[14];
			breathingEst_xCorr_mean[array_index] = (double)vsos_array[15];
			breath_ti[array_index] = (double)vsos_array[25];
			heart_ti[array_index] = (double)vsos_array[26];
		}

        // current_window_bmi
        if (array_index_bmi < 1200) {
            current_window_bmi[array_index_bmi] = (double)vsos_array[7];
			array_index_bmi++;
        }
        else {
            for (int num = 0; num < 1199; num++)
                current_window_bmi[num] = current_window_bmi[num + 1];
            current_window_bmi[1199] = (double)vsos_array[7];
        }

        // LF_HF_LFHF_windows
        if (array_index < 6000) {
            LF_HF_LFHF_windows[array_index] = (double)vsos_array[7];
            array_index++;
        }
        else {
            for (int num = 0; num < 5999; num++)
                LF_HF_LFHF_windows[num] = LF_HF_LFHF_windows[num + 1];
            LF_HF_LFHF_windows[5999] = (double)vsos_array[7];
        }

        // Main
		if (array_index >= 800) {
			// Array shift
			array_shift();
			unwrapPhasePeak_mm[799] = (double)vsos_array[7];
			heartRateEst_FFT_mean[799] = (double)vsos_array[10];
			heartRateEst_xCorr_mean[799] = (double)vsos_array[12];
			breathingEst_FFT_mean[799] = (double)vsos_array[14];
			breathingEst_xCorr_mean[799] = (double)vsos_array[15];
			breath_ti[799] = (double)vsos_array[25];
			heart_ti[799] = (double)vsos_array[26];

			// Setting time (initial)
			time(&now);
			struct tm *local = localtime(&now);
			hours = local->tm_hour;         // 獲取自午夜以來的小時數 (0-23)
			minutes = local->tm_min;        // 獲取小時後經過的分鐘數 (0-59)
			seconds = local->tm_sec;        // 獲取一分鐘後經過的秒數 (0-59)
			day = local->tm_mday;            // 獲取月份中的日期(1 到 31)
			month = local->tm_mon + 1;      // 獲取一年中的月份(0 到 11)
			year = local->tm_year + 1900;   // 獲取自 1900 年以來的年份
			
			// Setting time (define)
			if (year - start_year >= 1){
				start_year = year;
				start_month = month;
				start_day = day;
				start_time_sec = seconds+minutes*60+hours*3600;
				next_YMD = 1;
			}
			else if (month - start_month >= 1){
				start_month = month;
				start_day = day;
				start_time_sec = seconds+minutes*60+hours*3600;
				next_YMD = 1;
			}
			else if (day - start_day >= 1){
				start_day = day;
				start_time_sec = seconds+minutes*60+hours*3600;
				next_YMD = 1;
			}

			if (end_time - start_time >= 1 || next_YMD == 1)
			{
                
				next_YMD = 0;

                // Reset average parameters
				hr_mean_FFT = 0;
				hr_mean_xCorr = 0;
				br_mean_FFT = 0;
				br_mean_xCorr = 0;
				breath_mean_ti = 0;
				heart_mean_ti = 0;

				// Calculated average
				for (int i = 0; i < 800; i++){
					hr_mean_FFT += heartRateEst_FFT_mean[i];
					hr_mean_xCorr += heartRateEst_xCorr_mean[i];
					br_mean_FFT += breathingEst_FFT_mean[i];
					br_mean_xCorr += breathingEst_xCorr_mean[i];
					breath_mean_ti += breath_ti[i];
					heart_mean_ti += heart_ti[i];
				}
				hr_mean_FFT = hr_mean_FFT / 800;
				hr_mean_xCorr = hr_mean_xCorr / 800;
				br_mean_FFT = br_mean_FFT / 800;
				br_mean_xCorr = br_mean_xCorr / 800;
				breath_mean_ti = breath_mean_ti / 800;
				heart_mean_ti = heart_mean_ti / 800;

				/* ---------------------------- Breath Heart ---------------------------- */
				// --------------------- Phase_difference --------------------- 
				size_t len_unwrapPhasePeak_mm = sizeof(unwrapPhasePeak_mm) / sizeof(unwrapPhasePeak_mm[0]);
				for (int num = 1; num < len_unwrapPhasePeak_mm; num++)
				{
					phase_diff[num - 1] = unwrapPhasePeak_mm[num] - unwrapPhasePeak_mm[num - 1];
				}

				// For breathing heartbeat loop
				for (int br0hr1 = 0; br0hr1 < 2; br0hr1++){

					// Remove_impulse_noise
					size_t len_phase_diff = sizeof(phase_diff) / sizeof(phase_diff[0]);
					if (br0hr1 == 0)
						thr = 1.5;
					else
						thr = 1.5;
					for (int num = 1; num < len_phase_diff - 1; num++)
					{
						forward = phase_diff[num] - phase_diff[num - 1];
						backward = phase_diff[num] - phase_diff[num + 1];
						if ((forward > thr && backward > thr) || (forward < -thr && backward < -thr))
						{
							// printf("%f\n", phase_diff[num-1] + (phase_diff[num+1] -  phase_diff[num-1])/2);
							removed_noise[num] = (double)phase_diff[num - 1] + (double)(phase_diff[num + 1] - (double)phase_diff[num - 1]) / 2;
						}
						removed_noise[num] = (double)phase_diff[num];
					}

					// --------------------- iir_bandpass_filter_1 --------------------- 
					//  def zpk2tf(z, p, k):
					//  if (order == 5) => BR
					double y[799] = {0};
					if (br0hr1 == 0) {
						double b[11] = {0.0003100856139325839, -0.0024503932074526414, 0.008191253474064738, -0.014462792713750464, 0.012602969766102465, 0.0, -0.012602969766102468, 0.01446279271375046, -0.008191253474064738, 0.0024503932074526414, -0.00031008561393258384};
						double a[11] = {1., -9.780812442849507, 43.08317719449593, -112.54898060868825, 193.10308878043222, -227.3654181511635, 186.0550507759317, -104.48315327333019, 38.53588426508279, -8.42919504975534, 0.8303585098573365};
						double delay[10] = {0}; // length of (a or b) - 1
						size_t len_ab = sizeof(a) / sizeof(double);
						double after_lifter[799] = {0};
						lfilter(b, a, removed_noise, y, delay, len_ab, 799, 1, 1);
					}

					//  if (order == 9) => HR
					else {
						double b[19] = {0.0009309221423942934, -0.012651127859899214, 0.08140072903422219, -0.3279827818967048, 0.9206296690623369, -1.8881391347626006, 2.864155296830487, -3.1158825565644306, 2.0796297285426384, -1.6933375125414267e-15, -2.0796297285426384, 3.1158825565644315, -2.8641552968304875, 1.8881391347626006, -0.9206296690623373, 0.3279827818967048, -0.08140072903422219, 0.012651127859899216, -0.0009309221423942937};
						double a[19] = {1.0, -15.142612391789038, 109.52867216475022, -502.6626423702458, 1639.7914526180646, -4037.095860679349, 7772.407643496967, -11963.349399532222, 14922.878739349395, -15197.385490007382, 12665.576290591296, -8617.837327643654, 4751.9970631301, -2094.9246186588152, 722.2258948681153, -187.91289244221497, 34.75519411187606, -4.078772963228582, 0.22866640874033417};
						double delay[18] = {0};  // length of (a or b) - 1
						size_t len_ab = sizeof(a) / sizeof(double);
						double after_lifter[799] = {0};
						lfilter(b, a, removed_noise, y, delay, len_ab, 799, 1, 1);
					}
					
					// --------------------- FFT --------------------- 
					int N = 799;
					double P[800];
					rfft_forward_1d_array(y, N, N, 1, 1, P);  // Output: y

					// Find max value index in FFT
					int max_index = 0;
					double max = sqrt(pow(P[0], 2) + pow(P[0 + 1], 2));

					for (int i = 0; i < 798; i += 2)
					{
						if (sqrt(pow(P[i], 2) + pow(P[i + 1], 2)) > max)
						{
							max = sqrt(pow(P[i], 2) + pow(P[i + 1], 2));
							max_index = i;
						}
					}
					double index_of_fftmax = (max_index / 2) * 10.0 / (int)(N / 2);  // Output: index_of_fftmax

					// --------------------- Smoothing signal --------------------- 
					int smoothing_pars;
					if (br0hr1 == 0)
						smoothing_pars = 2;
					else
						smoothing_pars = 2;
					int len_input = sizeof(y) / sizeof(double);
					double *data_s = MLR(y, smoothing_pars, len_input);  // Output: data_s
					
					// --------------------- Feature_detection ---------------------
					// Signal length and half length
					len_s_half = floor(len_input / 2);

					// Output peak
					m_p = 0;  // Pointer to the end of valid area in allocated arrays
					int midpoints_peak[len_s_half];
					feature_peak = _local_maxima_1d(data_s, len_input, midpoints_peak, &m_p);  // Output: feature_peak

					// Dynamic memory management for negative signals
					for (int i = 0; i < len_input; i++)
						neg_x[i] = -data_s[i];

					// Output valley
					m_v = 0;  // Pointer to the end of valid area in allocated arrays
					int midpoints_valley[len_s_half];
					feature_valley = _local_maxima_1d(neg_x, len_input, midpoints_valley, &m_v);  // Output: feature_valley
					
					// --------------------- Feature compress --------------------- 
					// Initialize (feature_compress)

					// Given value
					for (int i = 0; i < m_p; i++)
						total_feature[i] = feature_peak[i];
					for (int i = 0; i < m_v; i++)
						total_feature[m_p+i] = feature_valley[i];

					// Perform quicksort on data
					quickSort(total_feature, 0, m_p + m_v - 1);
					
					// Compress parameter
					if (br0hr1 == 0)
						time_thr = 22;
					else
						time_thr = 5;
					ltera = 0;
					sum_v = 0;
					sum_p = 0;

					// Algorithm
					while (ltera < m_p + m_v - 1)
					{
						// Record start at valley or peak (peak:0 valley:1)
						start_feature = 1;
						for (int i = 0; i < m_p; i++){
							if (feature_peak[i] == total_feature[ltera]){
								start_feature = 0;
								break;
							}
						}

						ltera_add = ltera;
						while (total_feature[ltera_add+1]-total_feature[ltera_add]<time_thr){
							// skip the feature which is too close
							ltera_add += 1;
							// break the loop if it is out of boundary
							if(ltera_add >= (m_p + m_v - 1))
								break;
						}

						// Record end at valley or peak (peak:0 valley:1)
						end_feature = 1;
						for (int i = 0; i < m_p; i++){
							if (feature_peak[i] == total_feature[ltera_add]){
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
					printf("\n");
					
					// --------------------- Feature sort --------------------- 

					// Given value
					for (int i = 0; i < sum_p; i++)
						compress_feature[i] = feature_compress_peak[i];
					for (int i = 0; i < sum_v; i++)
						compress_feature[sum_p+i] = feature_compress_valley[i];
					
					// Perform quicksort on data
					quickSort(compress_feature, 0, sum_p + sum_v - 1);
					
					// --------------------- Candidate search --------------------- 
					double tmp_sum, window_sum, tmp_var, window_var;
					int NT_index, NB_index;

					// Doing the zero paddding
					int window_size;
					if (br0hr1 == 0)
						window_size = 17;
					else
						window_size = 4;

					for (int i = 0; i < len_input + (2 * window_size); i++){
						if (i >= window_size && i < len_input + (2 * window_size) - window_size)
							signal_pad[i] = data_s[i-window_size];
						else
							signal_pad[i] = 1;
					}
					for (int i = 0; i < window_size; i++)
						signal_pad[i] = data_s[0];
					for (int i = len_input + (2 * window_size) - window_size; i < len_input + (2 * window_size) - 1; i++)
						signal_pad[i] = data_s[len_input-1];

					// Calaulate the mean and std using windows(for peaks)
					NB_index = 0;
					NT_index = 0;
					for (int i = 0; i < sum_v+sum_p; i++){
						// For the mean
						tmp_sum = 0;
						window_sum = 0;
						tmp_var = 0;
						window_var = 0;
						for (int j = compress_feature[i]; j < compress_feature[i]+2*window_size+1; j++)
							tmp_sum += signal_pad[j];
						window_sum = tmp_sum / (window_size*2+1);
						for (int j = compress_feature[i]; j < compress_feature[i]+2*window_size+1; j++)
							tmp_var += pow(signal_pad[j] - window_sum, 2) ;
						window_var = sqrt(tmp_var / (window_size*2+1));

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
					// --------------------- Caculate breath rate --------------------- 
					double rate, cur_rate, tmp_rate;
					rate = 0;
					cur_rate = 0;
					tmp_rate = 0;
					// If both NT and NB are not detected
					if (NT_index <= 1 && NB_index <= 1)
						rate = 0;

					// If only NT are detected
					else if (NT_index > 1 && NB_index <= 1){
						for (int i = 1; i < NT_index; i++)
							tmp_rate += NT_point[i] - NT_point[i-1];
						rate = 1200 / (tmp_rate / (NT_index - 1));
					}

					// If only NB are detected
					else if (NT_index <= 1 && NB_index > 1){
						for (int i = 1; i < NB_index; i++)
							tmp_rate += NB_point[i] - NB_point[i-1];
						rate = 1200 / (tmp_rate / (NB_index - 1));
					}

					// If both NT and NB are detected
					else {
						for (int i = 1; i < NT_index; i++)
							tmp_rate += NT_point[i] - NT_point[i-1];
						cur_rate = tmp_rate / (NT_index - 1);
						tmp_rate = 0;
						for (int i = 1; i < NB_index; i++)
							tmp_rate += NB_point[i] - NB_point[i-1];
						cur_rate += tmp_rate / (NB_index - 1);
						rate = 1200 / (cur_rate / 2);
					}

					// SVC model
					svm_input[0] = index_of_fftmax;
					if (br0hr1 == 0){
						svm_input[1] = br_mean_FFT;
						svm_input[2] = br_mean_xCorr;
						svm_result = predict_br(svm_input);
						if (svm_result == 0)
							final_br = rate;
						else
							final_br = breath_mean_ti;
					}
					else{
						svm_input[1] = br_mean_FFT;
						svm_input[2] = br_mean_xCorr;
						svm_result = predict_hr(svm_input);
						if (svm_result == 0)
							final_hr = rate;
						else
							final_hr = heart_mean_ti;
					}
				}

				// 當呼吸律或心律為異常值, 以前一秒的輸出值取代。
				br_rpm = final_br;
				hr_rpm = final_hr;
				br_rpm = substitute(tmp_br, br_rpm, 1);
				hr_rpm = substitute(tmp_br, hr_rpm, 0);
				br_rpm = round(br_rpm*10000)/10000;
				hr_rpm = round(hr_rpm*10000)/10000;
				tmp_br = br_rpm;
				tmp_hr = hr_rpm;

				printf("BR = %f\nHR = %f\n", br_rpm, hr_rpm);
				
				if (seconds == 0 && counter == 0) {
					counter += 1;
					begin = 1;
					printf("Start recording features\n");
				}

				if (begin == 1) {
					if (var_index >= 600) {
						for (int num = 0; num < 599; num++) {
							var_RPM_br_KNN[num] = var_RPM_br_KNN[num + 1];
							var_RPM_hr_KNN[num] = var_RPM_hr_KNN[num + 1];
						}
						var_RPM_br_KNN[599] = br_rpm;
						var_RPM_hr_KNN[599] = hr_rpm;
					}
					else {
						var_RPM_br_KNN[var_index] = br_rpm;
						var_RPM_hr_KNN[var_index] = hr_rpm;
						var_index += 1;
					}

					// mov_dens (window size = 60 * 20)
					if (array_index_bmi == 1200) {

						// Function
						mov_dens_fn(current_window_bmi, &mov_dens);

						// Output (with slide)
						if (mov_dens_index >= 60) {
							for (int num = 0; num < 59; num++)
								mov_dens_ar[num] = mov_dens_ar[num + 1];
							mov_dens_ar[59] = mov_dens;
						}
						else {
							mov_dens_ar[mov_dens_index] = mov_dens;
							mov_dens_index++;
						}
					}

					/* ----------------------- Breath ----------------------- */
					// tfRSA
					if (var_index >= 10) {

						// Input
						for (int num = 9; num >= 0; num--)
							tfRSA_contain10_br[9-num] = var_RPM_br_KNN[var_index-1-num];
						
						// Function
						tfRSA_fn(tfRSA_contain10_br, &tfRSA);

						// Output (with slide)
						if (tfRSA_index >= 31) {
							for (int num = 0; num < 30; num++)
								tfRSA_ar[num] = tfRSA_ar[num + 1];
							tfRSA_ar[30] = tfRSA;
							open_stfRSA = 1;
						}
						else {
							tfRSA_ar[tfRSA_index] = tfRSA;
							tfRSA_index++;
						}
					}

					// sfRSA & sdfRSA
					if (var_index >= 31) {

						// ----------- Input (sfRSA) ----------- 
						for (int num = 30; num >= 0; num--)
							sfRSA_contain31_br[30-num] = var_RPM_br_KNN[var_index-1-num];
						
						// Function
						savgol_filter(sfRSA_contain31_br, 31, 3, &sfRSA_mean, sfRSA);

						// Output (with slide)
						if (sfRSA_index >= 60) {
							for (int num = 0; num < 59; num++)
								sfRSA_ar[num] = sfRSA_ar[num + 1];
							sfRSA_ar[59] = sfRSA_mean;
						}
						else {
							sfRSA_ar[sfRSA_index] = sfRSA_mean;
							sfRSA_index++;
						}

						// ----------- Input (sdfRSA) ----------- 
						for (int tmp_idx = 0; tmp_idx < 31; tmp_idx++)
							sdfRSA_tmp[tmp_idx] = abs(var_RPM_br_KNN[tmp_idx] - sfRSA[tmp_idx]);

						// Function
						savgol_filter(sdfRSA_tmp, 31, 3, &sdfRSA_mean, sdfRSA);

						// Output (with slide)
						if (sdfRSA_index >= 60) {
							for (int num = 0; num < 59; num++)
								sdfRSA_ar[num] = sdfRSA_ar[num + 1];
							sdfRSA_ar[59] = sdfRSA_mean;
						}
						else {
							sdfRSA_ar[sdfRSA_index] = sdfRSA_mean;
							sdfRSA_index++;
						}
					}
					
					// stfRSA
					if (open_stfRSA == 1) {

						// Function
						savgol_filter(tfRSA_ar, 31, 2, &stfRSA_mean, stfRSA);

						// Output (with slide)
						if (stfRSA_index >= 60) {
							for (int num = 0; num < 59; num++)
								stfRSA_ar[num] = stfRSA_ar[num + 1];
							stfRSA_ar[59] = stfRSA_mean;
						}
						else {
							stfRSA_ar[stfRSA_index] = stfRSA_mean;
							stfRSA_index++;
						}
					}

					// NOT YET TESTING
					/* ----------------------- Heart ----------------------- */
					// tmHR
					if (var_index >= 10) {

						// Input
						for (int num = 9; num >= 0; num--)
							tmHR_contain10_hr[9-num] = var_RPM_hr_KNN[var_index-1-num];

						// Function
						tfRSA_fn(tmHR_contain10_hr, &tmHR);
						
						// Output (with slide)
						if (tmHR_index >= 31) {
							for (int num = 0; num < 30; num++)
								tmHR_ar[num] = tmHR_ar[num + 1];
							tmHR_ar[30] = tmHR;
							open_stmHR = 1;
						}
						else {
							tmHR_ar[tmHR_index] = tmHR;
							tmHR_index++;
						}
					}

					// smHR & sdmHR
					if (var_index >= 31) {

						// ----------- Input (smHR) ----------- 
						for (int num = 30; num >= 0; num--)
							smHR_contain31_hr[30-num] = var_RPM_hr_KNN[var_index-1-num];
						
						// Function
						savgol_filter(smHR_contain31_hr, 31, 3, &smHR_mean, smHR);

						// Output (with slide)
						if (smHR_index >= 60) {
							for (int num = 0; num < 59; num++)
								smHR_ar[num] = smHR_ar[num + 1];
							smHR_ar[59] = smHR_mean;
						}
						else {
							smHR_ar[smHR_index] = smHR_mean;
							smHR_index++;
						}

						// ----------- Input (sdmHR) ----------- 
						for (int tmp_idx = 0; tmp_idx < 31; tmp_idx++)
							sdmHR_tmp[tmp_idx] = abs(var_RPM_hr_KNN[tmp_idx] - smHR[tmp_idx]);

						// Function
						savgol_filter(sdmHR_tmp, 31, 3, &sdmHR_mean, sdmHR);

						// Output (with slide)
						if (sdmHR_index >= 60) {
							for (int num = 0; num < 59; num++)
								sdmHR_ar[num] = sdmHR_ar[num + 1];
							sdmHR_ar[59] = sdmHR_mean;
						}
						else {
							sdmHR_ar[sdmHR_index] = sdmHR_mean;
							sdmHR_index++;
						}
					}

					// stmHR
					if (open_stmHR == 1) {

						// Function
						savgol_filter(stmHR_ar, 31, 2, &stmHR_mean, stmHR);

						// Output (with slide)
						if (stmHR_index >= 60) {
							for (int num = 0; num < 59; num++)
								stmHR_ar[num] = stmHR_ar[num + 1];
							stmHR_ar[59] = stmHR_mean;
						}
						else {
							stmHR_ar[stmHR_index] = stmHR_mean;
							stmHR_index++;
						}
					}
					
					// LF_HF_LFHF
					emerge_sum_LF = 0;
					emerge_sum_HF = 0;
					for (int state_HL = 0; state_HL < 2; state_HL++) {
						if (state_HL == 0){
							double b[5] = { 0.0009995048070486716, -0.003994449797763787, 0.005989890331755259, -0.003994449797763788, 0.0009995048070486716 };
							double a[5] = { 1.0, -3.996631500101437, 5.9910822943901785, -3.9922680954261405, 0.9978176514624266 };
							double delay[4] = {0}; // length of (a or b) - 1
							size_t len_ab = sizeof(a) / sizeof(double);
							lfilter(b, a, LF_HF_LFHF_windows, input_energe, delay, len_ab, array_index, 1, 1);
							
							// energe
							rfft_forward_1d_array(input_energe, array_index, array_index, 1, 1, out_fft);
							for (int i = 0; i < array_index-1; i+=2){
								emerge_sum_LF += pow(out_fft[i], 2) + pow(out_fft[i+1], 2);
							}
						}
						else if (state_HL == 1){
							double b[5] = { 0.0010005991658251364, -0.003978263207457302, 0.005955363064995437, -0.003978263207457302, 0.0010005991658251364 };
							double a[5] = { 1.0, -3.9832035773796473, 5.961516438161199, -3.973322837534953, 0.995044958484506 };
							double delay[4] = {0}; // length of (a or b) - 1
							size_t len_ab = sizeof(a) / sizeof(double);
							lfilter(b, a, LF_HF_LFHF_windows, input_energe, delay, len_ab, array_index, 1, 1);

							// energe
							rfft_forward_1d_array(input_energe, array_index, array_index, 1, 1, out_fft);
							for (int i = 0; i < array_index-1; i+=2){
								emerge_sum_HF += pow(out_fft[i], 2) + pow(out_fft[i+1], 2);
							}
						}
					}
					LFHF_eng = emerge_sum_HF / emerge_sum_LF;

					// Output (with slide) [LF, HF, LFHF]
					if (LF_HF_LFHF_index >= 60) {
						for (int num = 0; num < 59; num++) {
							LF_ar[num] = LF_ar[num + 1];
							HF_ar[num] = HF_ar[num + 1];
							LFHF_ar[num] = LFHF_ar[num + 1];
						}
						LF_ar[59] = emerge_sum_LF;
						HF_ar[59] = emerge_sum_HF;
						LFHF_ar[59] = LFHF_eng;
					}
					else {
						LF_ar[LF_HF_LFHF_index] = emerge_sum_LF;
						HF_ar[LF_HF_LFHF_index] = emerge_sum_HF;
						LFHF_ar[LF_HF_LFHF_index] = LFHF_eng;
						LF_HF_LFHF_index++;
					}

					// Output (with slide) [HF_arr, LFHF_arr]
					if (sF_index >= 31) {
						for (int num = 0; num < 30; num++) {
							HF_arr[num] = HF_arr[num + 1];
							LFHF_arr[num] = LFHF_arr[num + 1];
						}
						HF_arr[30] = emerge_sum_HF;
						LFHF_arr[30] = LFHF_eng;
						open_HF = 1;
					}
					else {
						HF_arr[sF_index] = emerge_sum_HF;
						LFHF_arr[sF_index] = LFHF_eng;
						sF_index++;
					}

					if (open_HF == 1) {

						// Function
						// sHF
						savgol_filter(HF_arr, 31, 3, &sHF_mean, sHF);
						// sLFHF
						savgol_filter(LFHF_arr, 31, 3, &sLFHF_mean, sLFHF);

						// Output (with slide)
						if (sHF_index >= 60) {
							for (int num = 0; num < 59; num++) {
								sHF_ar[num] = sHF_ar[num + 1];
								sLFHF_ar[num] = sLFHF_ar[num + 1];
							}
							sHF_ar[59] = sHF_mean;
							sLFHF_ar[59] = sLFHF_mean;
						}
						else {
							sHF_ar[sHF_index] = sHF_mean;
							sLFHF_ar[sHF_index] = sLFHF_mean;
							sHF_index++;
						}
					}

					// ------------------------------------- 睡眠階段 (Paper2) ------------------------------------- 
					// Variance of RPM
					if (var_index == 600) {
						var_RPM(var_RPM_br_KNN, &var_RPM_br);
						var_RPM(var_RPM_hr_KNN, &var_RPM_hr);

						// Output (with slide) 
						if (var_RPM_index >= 60) {
							for (int num = 0; num < 59; num++) {
								var_RPM_ar[num] = var_RPM_ar[num + 1];
								var_HPM_ar[num] = var_HPM_ar[num + 1];
							}
							var_RPM_ar[59] = var_RPM_br;
							var_HPM_ar[59] = var_RPM_hr;
						}
						else {
							var_RPM_ar[var_RPM_index] = var_RPM_br;
							var_HPM_ar[var_RPM_index] = var_RPM_hr;
							var_RPM_index++;
						}
						looper += 1;
					}

					// Body Movement Index (BMI)
                	// Deep Parameter
					if (array_index_bmi == 1200) {
						bmi(current_window_bmi, &bmi_current);
						hk = 0;
						for (int i = var_index-60; i < var_index; i ++)
							hk += var_RPM_hr_KNN[i];
						hk = hk / 60;
						deep_parameter(bmi_current, hk, &dk);

						// Output (with slide) 
						if (bmi_ar_index >= 60) {
							for (int num = 0; num < 59; num++) {
								bmi_ar[num] = bmi_ar[num + 1];
								deep_p_ar[num] = deep_p_ar[num + 1];
							}
							bmi_ar[59] = bmi_current;
							deep_p_ar[59] = dk;
						}
						else {
							bmi_ar[bmi_ar_index] = bmi_current;
							deep_p_ar[bmi_ar_index] = dk;
							bmi_ar_index++;
						}

						// Amplitude Difference Accumulation (ADA) of Respiration
						double output_br = 0;
						double output_hr = 0;
    					ada(current_window_bmi, array_index_bmi, 0, &output_br);
						ada(current_window_bmi, array_index_bmi, 1, &output_hr);

						// Output (with slide) 
						if (ada_index >= 60) {
							for (int num = 0; num < 59; num++) {
								ada_br_ar[num] = ada_br_ar[num + 1];
								ada_hr_ar[num] = ada_hr_ar[num + 1];
							}
							ada_br_ar[59] = output_br;
							ada_hr_ar[59] = output_hr;
						}
						else {
							ada_br_ar[ada_index] = output_br;
							ada_hr_ar[ada_index] = output_hr;
							ada_index++;
						}
					}

					// REM Parameter
					if (var_index >= 5 * 60) {
						for (int i = 0; i < 300; i++)
							rem_contain[i] = var_RPM_br_KNN[var_index - (299 - i)];
						rem_parameter(rem_contain, &rem_par);

						// Output (with slide) 
						if (rem_par_index >= 60) {
							for (int num = 0; num < 59; num++) {
								rem_parameter_ar[num] = rem_parameter_ar[num + 1];
							}
							rem_parameter_ar[59] = rem_par;
						}
						else {
							rem_parameter_ar[rem_par_index] = rem_par;
							rem_par_index++;
						}
					}

					// Time features
					time(&now_tf);
					struct tm *local_tf = localtime(&now_tf);
					hours_tf = local_tf->tm_hour;
					minutes_tf = local_tf->tm_min;
					seconds_tf = local_tf->tm_sec;
					time_fn(hours_tf, minutes_tf, seconds_tf, &time_featres);

					// Output (with slide) 
					if (time_fn_index >= 60) {
						for (int num = 0; num < 59; num++) {
							time_ar[num] = time_ar[num + 1];
						}
						time_ar[59] = time_featres;
					}
					else {
						time_ar[time_fn_index] = time_featres;
						time_fn_index++;
					}
					
					// breath & heart
					// Output (with slide) 
					if (br_hr_index >= 60) {
						for (int num = 0; num < 59; num++) {
							breath_ar[num] = breath_ar[num + 1];
							heart_ar[num] = heart_ar[num + 1];
						}
						heart_ar[59] = hr_rpm;
						breath_ar[59] = br_rpm;
					}
					else {
						heart_ar[br_hr_index] = hr_rpm;
						breath_ar[br_hr_index] = br_rpm;
						br_hr_index++;
					}

					// 換日
					if (counter_mean == 0) {
						start_hour = hours_tf;
						start_min = minutes_tf;
						end_hour = start_hour;
						end_min = start_min;
						counter_mean += 1;
					}
					else {
						end_hour = hours_tf;
						end_min = minutes_tf;
					}

					// 小時
					if (end_hour - start_hour >= 1 || end_hour - start_hour == -23) {
						start_hour = end_hour;
						start_min = end_min;
						next_HM = 1;
					}

					// 分鐘
					if (end_min - start_min >= 1) {
						start_min = end_min;
						next_HM = 1;
					}

					if (next_HM == 1 && looper >= 10) {
						mean_fn_double (breath_ar, br_hr_index, &all_results[0]);
						mean_fn_double (heart_ar, br_hr_index, &all_results[1]);
						mean_fn_double (bmi_ar, bmi_ar_index, &all_results[2]);
						mean_fn_double (deep_p_ar, bmi_ar_index, &all_results[3]);
						mean_fn_double (ada_br_ar, ada_index, &all_results[4]);
						mean_fn_double (ada_hr_ar, ada_index, &all_results[5]);
						mean_fn_double (var_RPM_ar, var_RPM_index, &all_results[6]);
						mean_fn_double (var_HPM_ar, var_RPM_index, &all_results[7]);
						mean_fn_double (rem_parameter_ar, rem_par_index, &all_results[8]);
						mean_fn_double (mov_dens_ar, mov_dens_index, &all_results[9]);
						mean_fn_double (LF_ar, LF_HF_LFHF_index, &all_results[10]);
						mean_fn_double (HF_ar, LF_HF_LFHF_index, &all_results[11]);
						mean_fn_double (LFHF_ar, LF_HF_LFHF_index, &all_results[12]);
						mean_fn_double (sHF_ar, sHF_index, &all_results[13]);
						mean_fn_double (sLFHF_ar, sHF_index, &all_results[14]);
						mean_fn_double (tfRSA_ar, tfRSA_index, &all_results[15]);
						mean_fn_double (tmHR_ar, tmHR_index, &all_results[16]);
						mean_fn_double (sfRSA_ar, sfRSA_index, &all_results[17]);
						mean_fn_double (smHR_ar, smHR_index, &all_results[18]);
						mean_fn_double (sdfRSA_ar, sdfRSA_index, &all_results[19]);
						mean_fn_double (sdmHR_ar, sdmHR_index, &all_results[20]);
						mean_fn_double (stfRSA_ar, stfRSA_index, &all_results[21]);
						mean_fn_double (stmHR_ar, stmHR_index, &all_results[22]);
						mean_fn_int (time_ar, time_fn_index, &all_results[23]);

						for (int re_index = 0; re_index < 24; re_index++) {
							if (all_results[re_index] > 140700000)
							    all_results[re_index] = 140700000;
						}

						// Random forest classifier prediction
						predict_result = predict(all_results);

						// Write data to csv
						FILE *fp = fopen(filename, "a");
						if (fp == NULL)
						{
							printf("error");
							return -1;
						}
						fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d:%d:%d, %d\n", all_results[1], all_results[0], all_results[2], all_results[3], all_results[4], all_results[5], all_results[6], all_results[7], all_results[8], all_results[9], all_results[10], all_results[11], all_results[12], all_results[13], all_results[14], all_results[15], all_results[16], all_results[17], all_results[18], all_results[19], all_results[20], all_results[21], all_results[22], all_results[23], hours_tf, minutes_tf, seconds_tf, predict_result);
						fclose(fp);

						// Parameter initialization
						br_hr_index = 0;
						bmi_ar_index = 0;
						ada_index = 0;
						var_RPM_index = 0;
						rem_par_index = 0;
						mov_dens_index = 0;
						LF_HF_LFHF_index = 0;
						sHF_index = 0;
						tfRSA_index = 0;
						tmHR_index = 0;
						sfRSA_index = 0;
						smHR_index = 0;
						sdfRSA_index = 0;
						sdmHR_index = 0;
						stfRSA_index = 0;
						stmHR_index = 0;
						time_fn_index = 0;
						next_HM = 0;

					}
				}

				start_time = end_time;
				printf("\n");
            }
        }
    }
    close(serial_port);
    return 0;
}