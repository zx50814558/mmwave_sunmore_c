// sudo chmod 666 ttyTHS1
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
#include "pocketfft.h"

// sklearn model
#include "svm_br_office_all.h"
#include "svm_hr_office_all.h"

struct timeval start, stop;
double secs = 0;

/*
input: 輸入訊號
delta: 左右延伸的窗格長度
len_: 輸入資料長度*/
double *MLR(double *input, int delta, int len_)
{
	double *data_s = input;  // 取得輸入訊號

	// 建立空矩陣
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

			// 窗格內算平均
			double mean_tmp = 0;
			for (int j = 0; j < end - star; j++)
				mean_tmp += input[star + j];
			mean_ptr[t] = (mean_tmp / (end - star));

			// Slope & bias
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
				tmp_s += (m_ptr[i] * t) + b_ptr[i];  // y = mt + b
			}
			data_s[t] = tmp_s / (2 * delta + 1);  // mean
		}
	}

	return data_s;
}

/*
x: 輸入訊號
len_s: 輸入訊號長度，因為是平滑化後的值，在尾端加 s
midpoints: 以指標傳入的陣列，用來回傳值
m: 指標，用以記錄 local maxima features 的數量*/
int *_local_maxima_1d(double *x, int len_s, int *midpoints, int *m)
{

	// 初始化陣列
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
		i += 1;  // Next elements
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

/*
array: 輸入陣列
low: Sort 最小位置 ( 開始 )
high: Sort 最小位置 ( 結束 )*/
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

// Data containers
float unwrapPhasePeak_mm[800];
float heartRateEst_FFT_mean[800];
float heartRateEst_xCorr_mean[800];
float breathingEst_FFT_mean[800];
float breathingEst_xCorr_mean[800];
float breath_ti[800];
float heart_ti[800];

// Sklearn to c
int svm_result;
double svm_input[3];
// double tmp_breath_rate = 0;

// Shifting data
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

/*Filter data along one-dimension with an IIR or FIR filter.
b: The numerator coefficient vector in a 1-D sequence.
a: The denominator coefficient vector in a 1-D sequence.
x: Signal after noise removal
y: Output filtered signal.
Z: Initial conditions for the filter delays. It is a vector (or array of vectors for an N-dimensional input) of length max(len(a), len(b)) - 1.
len_b: The parameter length of b and a.
len_x: The length of x
stride_X: default 1
stride_Y: default 1*/
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

void mean_fn_double (double *sig, int sig_len, double *output_d)
{
    *output_d = 0;
    for (int i = 0; i < sig_len; i++)
        *output_d += sig[i];
    *output_d = *output_d / sig_len;
}

int main(void)
{
	/* Initialize (Feature_detection) */
	int len_s_half;  // Signal length and half length
	int *feature_peak, *feature_valley;  // 回傳 Feature_detection 後 peak array 與 valley array 的指標
	int m_p, m_v; // Number of peak & valley
	double neg_x[799] = {0};  // 將原先的訊號上下反轉，用來取 feature_valley 的數值

	/* Initialize (feature_compress) */
	int start_feature, end_feature;  // Record start and end at valley or peak (peak:0 valley:1)
	int ltera_add;  // 用來做窗格內值的 Shift
	int time_thr;  // 判斷訊號間的 valley or peak 是否過於接近 (呼吸: 22, 心跳: 5)
	int location;  // 用於取值
	int ltera, sum_v, sum_p;  // 第幾輪, compress 後共幾個 valley, compress 後共幾個 peak

	int serial_port = open("/dev/ttyTHS1", O_RDWR);  // 設定 port 號
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

	/* Initialize */
	time_t start_time;  // 宣告時間變數 (開始時間)
	start_time = time(NULL);  // 讀取當前時間做為 (開始時間)
	int array_index = 0;  // 輸入值累加數量 ( 需累加到 800 個值才開始執行後續算法)
	float phase_diff[799];  // 宣告存放相位差的陣列 (長度 799)
	double removed_noise[799];  // 宣告存放 Remove_impulse_noise 後的陣列 (長度 799)
	float forward, backward;  // 算法中暫存的變數 ( forward: 窗格中心點減其一個值, backward: 窗格中心點減其一個值)
	float hr_mean_FFT, hr_mean_xCorr, br_mean_FFT, br_mean_xCorr, breath_mean_ti, heart_mean_ti;  // 用於計算平均數值，每項共 800 個值，累加後取平均
	float thr;  // 宣告 Remove_impulse_noise 時要濾除的閾值

	/* Feature compress */
	int total_feature[800] = {0};  // 存放所有的特徵 ( valley or peak )
	int feature_compress_valley[800] = {0};  // 存放 compress 後留下來沒被濾除的 valley 特徵
	int feature_compress_peak[800] = {0};  // 存放 compress 後留下來沒被濾除的 peak 特徵
	int compress_feature[800] = {0};  // 用於代表 Feature compress 後輸出的所有特徵 ( valley or peak )
	int NT_point[800] = {0};  // 累加 compress 後剩下幾個 peak feature
	int NB_point[800] = {0};  // 累加 compress 後剩下幾個 valley feature

	/* Candidate search */
	double signal_pad[1600] = {0};  // 存放 Candidate search 算法過程中的值

    /* People detect */
    double current_window_ebr[60] = {0};  // 用於累加呼吸律的平均能量，判斷人在不再
    double current_window_ehr[60] = {0};  // 用於累加呼吸律的平均能量，判斷人在不再
    int eng_index = 0;  // 陣列的 index 方便記錄目前蒐集到幾個值 ( 呼吸心律共用 )

	/* Final results  */
	double br_rate, hr_rate;  // 紀錄呼吸律與心律
	int hours, minutes, seconds;  // 完成一輪後的當下時間，用於 log
	time_t now_record;  // 宣告時間變數

	/* File handling */
    char filename[100] = {0};  // 檔名
    char input_name[100];
    char *input_n;
    char *root_dir = "dataset/"; // 檔案路徑
    printf("Input file name = ");  // 輸入檔名的提示
    input_n = fgets(input_name, 100, stdin);  // 寫入檔名
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
	fprintf(fp, "Times, heart, breath\n");
	fclose(fp);

	/* Start execution of the algorithm */
	while (1)
	{
		/* 從雷達讀取檔案資料 */
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

		/* 將讀取到並解碼後的資料累加，並接續使用 */
		time_t end_time;  // 宣告結束時間
		end_time = time(NULL);  // 以當前時間當作結束時間

		/* 呼吸律與心律的能量蒐集，當超過 60 個時向左 Shift 並推疊最新的數值到陣列尾端 */
		if (eng_index < 60)
		{
			current_window_ebr[eng_index] = vsos_array[22];  // 小於 60 時以索引值讀入呼吸能量
			current_window_ehr[eng_index] = vsos_array[23];  // 小於 60 時以索引值讀入心律能量
			eng_index++;  // 索引值 +1 ，下一輪指到下一個位置
		}
		else {
			// 當推疊資料超過 60 個，先做左 Shift
			for (int i_eng = 0; i_eng < 59; i_eng++) {
				current_window_ebr[i_eng] = current_window_ebr[i_eng + 1];  // 呼吸
				current_window_ehr[i_eng] = current_window_ehr[i_eng + 1];  // 心律
			}
			// 將新的值取代最後一個值
			current_window_ebr[59] = vsos_array[22];  // 呼吸
			current_window_ehr[59] = vsos_array[23];  // 心律
		}

		/* 演算法所需資料蒐集，當超過 800 個時向左 Shift 並推疊最新的數值到陣列尾端 */
		if (array_index < 800)
		{
			unwrapPhasePeak_mm[array_index] = vsos_array[7];  // unwrapPhasePeak_mm
			heartRateEst_FFT_mean[array_index] = vsos_array[10];  // heartRateEst_FFT
			heartRateEst_xCorr_mean[array_index] = vsos_array[12];  // heartRateEst_xCorr
			breathingEst_FFT_mean[array_index] = vsos_array[14];  // breathingEst_FFT
			breathingEst_xCorr_mean[array_index] = vsos_array[15];  // breathingEst_xCorr
			breath_ti[array_index] = vsos_array[25];  // ti 預測呼吸律
			heart_ti[array_index] = vsos_array[26];  // ti 預測心律
			array_index++;  // 更新索引值
		}
		else
		{
			array_shift();  // 左 Shift 矩陣內部所有值，進入這部分代表已經取得所需特徵數量，因此需更新陣列內所有值
			unwrapPhasePeak_mm[799] = vsos_array[7];  // 並推疊最新的數值到陣列尾端
			heartRateEst_FFT_mean[799] = vsos_array[10];  // 並推疊最新的數值到陣列尾端
			heartRateEst_xCorr_mean[799] = vsos_array[12];  // 並推疊最新的數值到陣列尾端
			breathingEst_FFT_mean[799] = vsos_array[14];  // 並推疊最新的數值到陣列尾端
			breathingEst_xCorr_mean[799] = vsos_array[15];  // 並推疊最新的數值到陣列尾端
			breath_ti[799] = vsos_array[25];  // 並推疊最新的數值到陣列尾端
			heart_ti[799] = vsos_array[26];  // 並推疊最新的數值到陣列尾端

			/* 當結束時間 - 開始時間 >= 1，開始執行算法，代表每間格 1 秒執行一次 */
			if (end_time - start_time >= 1)
			{

				/* 紀錄當前時間 */
				time(&now_record);
				struct tm *local = localtime(&now_record);
				hours = local->tm_hour;         // 獲取自午夜以來的小時數 (0-23)
				minutes = local->tm_min;        // 獲取小時後經過的分鐘數 (0-59)
				seconds = local->tm_sec;        // 獲取一分鐘後經過的秒數 (0-59)

				/* 計算參數平均 */
				// Reset average parameters
				hr_mean_FFT = 0;
				hr_mean_xCorr = 0;
				br_mean_FFT = 0;
				br_mean_xCorr = 0;
				breath_mean_ti = 0;
				heart_mean_ti = 0;
				// Calculated average ( 累加 )
				for (int i = 0; i < 800; i++){
					hr_mean_FFT += heartRateEst_FFT_mean[i];
					hr_mean_xCorr += heartRateEst_xCorr_mean[i];
					br_mean_FFT += breathingEst_FFT_mean[i];
					br_mean_xCorr += breathingEst_xCorr_mean[i];
					breath_mean_ti += breath_ti[i];
					heart_mean_ti += heart_ti[i];
				}
				// Calculated average ( 取平均 )
				hr_mean_FFT = hr_mean_FFT / 800;
				hr_mean_xCorr = hr_mean_xCorr / 800;
				br_mean_FFT = br_mean_FFT / 800;
				br_mean_xCorr = br_mean_xCorr / 800;
				breath_mean_ti = breath_mean_ti / 800;
				heart_mean_ti = heart_mean_ti / 800;

				// printf("\nhr_mean_FFT = %f\nhr_mean_xCorr = %f\nbr_mean_FFT = %f\nbr_mean_xCorr = %f\n", hr_mean_FFT, hr_mean_xCorr, br_mean_FFT, br_mean_xCorr);
				// printf("\nTI BR = %f, TI HR = %f\n", breath_mean_ti, heart_mean_ti);
				start_time = end_time;  // 當執行上述步驟後，更開始時間為結束時間，以便後續間隔 1 秒執行

				// --------------------- Phase_difference --------------------- 
				size_t len_unwrapPhasePeak_mm = sizeof(unwrapPhasePeak_mm) / sizeof(unwrapPhasePeak_mm[0]);  // 取得 unwrapPhasePeak_mm 長度
				for (int num = 1; num < len_unwrapPhasePeak_mm; num++)
				{
					phase_diff[num - 1] = unwrapPhasePeak_mm[num] - unwrapPhasePeak_mm[num - 1];  // 將後一個值減去前一個並存入 phase_diff 中
				}

				// For breathing heartbeat loop ( 0: 呼吸, 1: 心律)
				for (int br0hr1 = 0; br0hr1 < 2; br0hr1++){

					// Remove_impulse_noise
					size_t len_phase_diff = sizeof(phase_diff) / sizeof(phase_diff[0]);  // 取得 phase_diff 長度
					
					// 根據呼吸律或是心律設定不同的閾值
					if (br0hr1 == 0)
						thr = 1.5;
					else
						thr = 1.5;
					
					// 第一項無須處理，直接放入新的陣列
					removed_noise[0] = phase_diff[0];

					// removed_noise 演算法: 在窗格大小為 3 的陣列當中，中間項減去前一項 = forward, 中間項減去後一項 = backward
					for (int num = 1; num < len_phase_diff - 1; num++)
					{
						forward = phase_diff[num] - phase_diff[num - 1];  // 中間項減去前一項
						backward = phase_diff[num] - phase_diff[num + 1];  // 中間項減去後一項

						// 若兩項特徵同時超出閾值則設進行處理
						if ((forward > thr && backward > thr) || (forward < -thr && backward < -thr))
						{
							// printf("%f\n", phase_diff[num-1] + (phase_diff[num+1] -  phase_diff[num-1])/2);
							removed_noise[num] = (double)phase_diff[num - 1] + (double)(phase_diff[num + 1] - (double)phase_diff[num - 1]) / 2;  // 處理算法: 前一項 + (後一項 - 前一項) / 2
						}
						removed_noise[num] = (double)phase_diff[num];  // 若判定為正常值，沿用原始數據
					}

					// --------------------- iir_bandpass_filter_1 --------------------- 
					double y[799] = {0};  // 完成濾波後的輸出，訊號與輸入前等長
					//  def zpk2tf(z, p, k):
					//  if (order == 5) => BR
					// b & a: 對應呼吸心律的濾波參數設定
					// delay: Initial conditions for the filter delays. It is a vector (or array of vectors for an N-dimensional input) of length max(len(a), len(b)) - 1.
					if (br0hr1 == 0) {
						double b[11] = {0.000310085613932583790096353393,-0.002450393207452640498278384484,0.008191253474064732684190026646,-0.014462792713750451459309154245,0.012602969766102464777013381081,0.000000000000000000000000000000,-0.012602969766102459572842953150,0.014462792713750454928756106199,-0.008191253474064734418913502623,0.002450393207452640498278384484,-0.000310085613932583735886244769};
						double a[11] = {1.0,-9.780812442849507348796578298789,43.083177194495931416895473375916,-112.548980608688253823856939561665,193.103088780432216253757360391319,-227.365418151163510174228576943278,186.055050775931675843821722082794,-104.483153273330174215516308322549,38.535884265082792410339607158676,-8.429195049755342949993064394221,0.830358509857336501980284992896};
						double delay[10] = {0}; // length of (a or b) - 1
						size_t len_ab = sizeof(a) / sizeof(double);
						lfilter(b, a, removed_noise, y, delay, len_ab, 799, 1, 1);
					}

					//  if (order == 9) => HR
					else {
						double b[19] = {0.0009309221423942934, -0.012651127859899214, 0.08140072903422219, -0.3279827818967048, 0.9206296690623369, -1.8881391347626006, 2.864155296830487, -3.1158825565644306, 2.0796297285426384, -1.6933375125414267e-15, -2.0796297285426384, 3.1158825565644315, -2.8641552968304875, 1.8881391347626006, -0.9206296690623373, 0.3279827818967048, -0.08140072903422219, 0.012651127859899216, -0.0009309221423942937};
						double a[19] = {1.0, -15.142612391789038, 109.52867216475022, -502.6626423702458, 1639.7914526180646, -4037.095860679349, 7772.407643496967, -11963.349399532222, 14922.878739349395, -15197.385490007382, 12665.576290591296, -8617.837327643654, 4751.9970631301, -2094.9246186588152, 722.2258948681153, -187.91289244221497, 34.75519411187606, -4.078772963228582, 0.22866640874033417};
						double delay[18] = {0};  // length of (a or b) - 1
						size_t len_ab = sizeof(a) / sizeof(double);
						lfilter(b, a, removed_noise, y, delay, len_ab, 799, 1, 1);
					}
					
					// --------------------- FFT --------------------- 
					int N = 799;  // FFT length & The number of samples
					double P[800];  // Output signal(complex-value). The layout of elemens are: `nrows * ((fft_len / 2) + 1) * 2(real, img)
					rfft_forward_1d_array(y, N, N, 1, 1, P);  // Output: y

					// Find max value index in FFT
					int max_index = 0;  // 設定追蹤最大值的 index ，當發現新的最大值即更新
					double max = sqrt(pow(P[0], 2) + pow(P[0 + 1], 2));  // 計算複數最大值，預設第一項為最大後續更新

					// 遍歷所有頻域輸出，找出最大值為在的頻率，間隔設為 2 因為有實虛
					for (int i = 0; i < 798; i += 2)
					{
						// 若當前值大於最大值則進行取代
						if (sqrt(pow(P[i], 2) + pow(P[i + 1], 2)) > max)
						{
							max = sqrt(pow(P[i], 2) + pow(P[i + 1], 2));  // 計算複數最大值
							max_index = i;  // 更新追蹤最大值的 index
						}
					}
					double index_of_fftmax = (max_index / 2) * 10.0 / (int)(N / 2);  // Output: index_of_fftmax

					// --------------------- Smoothing signal --------------------- 
					int smoothing_pars;  // Smoothing signal 所需的參數，製作以當前值向左右延伸 smoothing_pars 形成的窗格
					if (br0hr1 == 0)
						smoothing_pars = 2;
					else
						smoothing_pars = 2;
					int len_input = sizeof(y) / sizeof(double);  // 計算輸入資料長度
					double *data_s = MLR(y, smoothing_pars, len_input);  // Output: data_s
					
					// --------------------- Feature_detection ---------------------
					// Signal length and half length
					len_s_half = floor(len_input / 2);

					// Output peak
					m_p = 0;  // Pointer to the end of valid area in allocated arrays
					int midpoints_peak[len_s_half];  // 初始化空間給 Function
					feature_peak = _local_maxima_1d(data_s, len_input, midpoints_peak, &m_p);  // Output: feature_peak ( 最大值集合 )

					// Reverse up-down array
					for (int i = 0; i < len_input; i++)
						neg_x[i] = -data_s[i];  // neg_x: 存放顛倒後的訊號，好以能夠重複使用 _local_maxima_1d 找出最小值

					// Output valley
					m_v = 0;  // Pointer to the end of valid area in allocated arrays
					int midpoints_valley[len_s_half];  // 初始化空間給 Function
					feature_valley = _local_maxima_1d(neg_x, len_input, midpoints_valley, &m_v);  // Output: feature_valley ( 最小值集合 )
					
					// --------------------- Feature compress --------------------- 
					// Given value
					for (int i = 0; i < m_p; i++)
						total_feature[i] = feature_peak[i];
					for (int i = 0; i < m_v; i++)
						total_feature[m_p+i] = feature_valley[i];

					// Perform quicksort on data
					quickSort(total_feature, 0, m_p + m_v - 1);
					
					// Compress parameter ( 呼吸: 22, 心律: 5 )
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

					// Given value ( 全部特徵放一起，後續對整個陣列做 Sort )
					for (int i = 0; i < sum_p; i++)
						compress_feature[i] = feature_compress_peak[i];
					for (int i = 0; i < sum_v; i++)
						compress_feature[sum_p+i] = feature_compress_valley[i];
					
					// Perform quicksort on data
					quickSort(compress_feature, 0, sum_p + sum_v - 1);
					
					// --------------------- Candidate search --------------------- 
					// Initialize (feature_compress)
					double tmp_sum, window_sum, tmp_var, window_var;
					int NT_index, NB_index;

					// 設定 Candidate search 所需的參數 ( window_size: 呼吸 = 17, 心律 = 4 )
					int window_size;
					if (br0hr1 == 0)
						window_size = 17;
					else
						window_size = 4;

					// Doing the zero paddding
					for (int i = 0; i < len_input + (2 * window_size); i++){
						if (i >= window_size && i < len_input + (2 * window_size) - window_size)
							signal_pad[i] = data_s[i-window_size];  // paddding 中間
						else
							signal_pad[i] = 1;  // 其餘 paddding 1
					}
					for (int i = 0; i < window_size; i++)
						signal_pad[i] = data_s[0];  // paddding 前面
					for (int i = len_input + (2 * window_size) - window_size; i < len_input + (2 * window_size) - 1; i++)
						signal_pad[i] = data_s[len_input-1];  // paddding 後面

					// Calaulate the mean and std using windows(for peaks)
					NB_index = 0;  // 作為索引值同時也作為計算多少 Bottom features 的計數器
					NT_index = 0;  // 作為索引值同時也作為計算多少 Top features 的計數器
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
							tmp_rate += NT_point[i] - NT_point[i-1];  // 特徵區間內的間隔總和
						rate = 1200 / (tmp_rate / (NT_index - 1));  // 特徵區間內的間隔總和除以區間個數 = 特徵平均距離 => 20(取樣頻率) * 60(秒) / 特徵平均距離
					}

					// If only NB are detected
					else if (NT_index <= 1 && NB_index > 1){
						for (int i = 1; i < NB_index; i++)
							tmp_rate += NB_point[i] - NB_point[i-1];  // 特徵區間內的間隔總和
						rate = 1200 / (tmp_rate / (NB_index - 1));  // 特徵區間內的間隔總和除以區間個數 = 特徵平均距離 => 20(取樣頻率) * 60(秒) / 特徵平均距離
					}

					// If both NT and NB are detected  ( 做法與上述相同但分開計算 )
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

					// The SVM classifier determines whether to use TI output or Ours algorithm output (0: Ours, 1: TI).
					svm_input[0] = index_of_fftmax;
					if (br0hr1 == 0){
						svm_input[1] = br_mean_FFT;
						svm_input[2] = br_mean_xCorr;
						svm_result = predict_br(svm_input);
						if (svm_result == 0)
							br_rate = round(rate);
						else
							br_rate = round(breath_mean_ti);
						
					}
					else{
						svm_input[1] = br_mean_FFT;
						svm_input[2] = br_mean_xCorr;
						svm_result = predict_hr(svm_input);
						if (svm_result == 0)
							hr_rate = round(rate);
						else
							hr_rate = round(heart_mean_ti);
					}
				}

				// 判斷雷達前是否有人存在 ( 若無人則以 0 取代呼吸律與心律 )
				double thr_br = 0;
				double thr_hr = 0;
				mean_fn_double(current_window_ebr, eng_index, &thr_br);
				mean_fn_double(current_window_ehr, eng_index, &thr_hr);
				if (thr_br <= 200 && thr_hr <= 30) {
					br_rate = 0;
					hr_rate = 0;
				}

				// 寫入 logs 檔案
				FILE *fp = fopen(filename, "a");
				if (fp == NULL)
				{
					printf("error");
					return -1;
				}
				// tmp_breath_rate = (ceil((int)hr_rate * 1.0 / 4) + (int)br_rate) / 2;
				fprintf(fp, "%d:%d:%d, %d, %d\n", hours, minutes, seconds, (int)hr_rate, (int)br_rate);
				fclose(fp);
				printf("\e[1;1H");
				int systemArb = system("clear");
				printf("=========================================================\n");
				printf("|      Version: V1.0                                    |\n");
				printf("=========================================================\n");
				printf("|      Heart rate: %d    |    Respiratory rate: %d      |\n", (int)hr_rate, (int)br_rate);
				printf("=========================================================\n");
			}
		}
	}
	close(serial_port);
	return 0; // success
}
