//sudo chmod 666 ttyTHS1
// C library headers
#include <stdio.h>
#include <string.h>
// Linux headers
#include <fcntl.h> // Contains file controls like O_RDWR
#include <errno.h> // Error integer and strerror() function
#include <termios.h> // Contains POSIX terminal control definitions
#include <unistd.h> // write(), read(), close()
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <complex.h>
#include "dbscan.c"
#include <stdlib.h>
int mean_count = 0;
int fall_lying_count = 0;
int compare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

int cmpfunc (const void * a, const void * b)
{
	return (*(float*)a - *(float*)b);
}

typedef union {
	float f;
	struct
	{
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} raw;
} myfloat;

unsigned int convertToInt(int* arr, int low, int high)
{
	unsigned f = 0, i;
	for (i = high; i >= low; i--) {
		f = f + arr[i] * pow(2, high - i);
	}
	return f;
}

unsigned int ieee754_convert(int* ieee)
{
	myfloat var;
	unsigned f = convertToInt(ieee, 9, 31);
	var.raw.mantissa = f;
	f = convertToInt(ieee, 1, 8);
	var.raw.exponent = f;
	var.raw.sign = ieee[0];
	return var.f;
}

struct tlvTypeInfo{
	int unitByte;
	int stateString;
	int sbyte;
	int dataByte;
	int retCnt;
	float nPoint;
	};
struct tlvTypeInfo changeit(struct tlvTypeInfo, int state, int count);
struct tlvTypeInfo changeit(struct tlvTypeInfo s, int state, int count)
{
	s.unitByte = 20;
	s.sbyte = 8;
	s.dataByte = 0;
	s.stateString = 3;
	if (state = 2)
	{
		s.unitByte = 20;
		s.sbyte = 8;
		s.dataByte = 8;
	}
	s.retCnt =  count - s.unitByte - s.sbyte;
	s.nPoint =  s.retCnt / s.sbyte;
	return (s);
}
int main() {
  
  time_t rawtime;
  struct tm *info;
  char buffer[80];
  time(&rawtime);
  info = localtime(&rawtime);
  printf("現在時間%s\n", asctime(info));
  char csv_name[10];
  printf("輸入檔名+.csv：");
  scanf("%s", csv_name);
  char *filename = csv_name;
  int serial_port = open("/dev/ttyTHS1", O_RDWR);
  struct termios tty;
  Struct output;
  //output = dbscan_output();
  if(tcgetattr(serial_port, &tty) != 0) {
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
  tty.c_iflag &= ~(IGNBRK|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|ICRNL); 
  tty.c_oflag &= ~OPOST;
  tty.c_oflag &= ~ONLCR;
  tty.c_cc[VTIME] = 10;   
  tty.c_cc[VMIN] = 0;
  
  cfsetispeed(&tty, B921600);
  cfsetospeed(&tty, B921600);
  
  if (tcsetattr(serial_port, TCSANOW, &tty) != 0) {
      printf("Error %i from tcsetattr: %s\n", errno, strerror(errno));
      return 1;
  }
  
  char read_buf [1024];
  int size;
  int magicWord[8] = {2, 1, 4, 3, 6, 5, 8, 7};
  unsigned long int header_reader_output[9];
  unsigned long int byte[4];
  unsigned int ieee[32];
  int k=0;
  
  while(1) 
  {
	  while ((size = read(serial_port, &read_buf, sizeof(read_buf)-1))>0)
	  {
		  struct tlvTypeInfo tlvTypeInfo_value;
		  struct tlvTypeInfo tlvTypeInfo_value1;
		  printf("----------start-----------\n");
		  int state = 0;
		  if (state == 0)
		  {
			  int same = 0;
			  for (int ix=0; ix< 8; ++ix) 
			  {		
					
					if (read_buf[ix] == magicWord[ix])
					{
						//printf("%d\n", read_buf[ix]);
						same++;
						if (same==8)
						{
							state = 1;
						}
					}
			  }		  
		  }
		  
		  if (state ==1)
		  {	
			  printf("----------header-----------\n");
			  for (int ix=8; ix< 48; ++ix) 
			  {
				  if (ix<44)
				  {
					  if (ix%4==0)
					  {
						  printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)) );
					  }			  
				  }
				  else
				  {
					  if (ix%2==0)
					  {	  
						  printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8)));
					  }  
					  if (ix==47)
					  {
						  state = 2;
					  }
				  }
			  }
		  }
		  
		  if (state ==2)
		  {
			  int tlvLength;
			  printf("----------TLV Header-----------\n");
			  for (int ix=48; ix< 56; ++ix) 
			  {
				  if (ix%4==0)
				  {
					  tlvLength = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)); 
					  //printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)) );
				  }
				  
				  if (ix==55)
				  {
					  tlvTypeInfo_value1 = changeit(tlvTypeInfo_value, state, tlvLength);
					  /*
					  printf("unitByte:%d\n", tlvTypeInfo_value1.unitByte);
					  printf("stateString:%d\n", tlvTypeInfo_value1.stateString);
					  printf("sbyte:%d\n", tlvTypeInfo_value1.sbyte);
					  printf("dataByte:%d\n", tlvTypeInfo_value1.dataByte);
					  printf("retCnt:%d\n", tlvTypeInfo_value1.retCnt);
					  printf("nPoint:%f\n", tlvTypeInfo_value1.nPoint);
					  */
					  state = tlvTypeInfo_value1.stateString;
				  }
				  
			  }
		  }
		  float V6_unit[5];
		  if (state == 3)
		  {
			  //5f
			  
			  for (int ix2=0; ix2<5; ++ix2) 
			  {
				  for (int ix1=59+(ix2*4); ix1> 55+(ix2*4); --ix1) 
				  {
					  //printf("%dix1:\n", ix1);
					  for (int ix=0; ix< 8; ++ix) 
					  {
						  ieee[k] = (read_buf[ix1] << ix & 0x80) >> 7;
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
				  V6_unit[ix2] = var.f;
				  printf("unit %f\n", var.f);				  
			  }	
			  state = 4;		  
		  }
		  
		  //printf("total point clude %f\n", tlvTypeInfo_value1.retCnt/8.0);
		  int total_point;
		  total_point = (int) tlvTypeInfo_value1.retCnt/8.0;
		  int v6_point_before[5];
		  float v6_point_after[5];
		  float v6_2d_output[total_point][5];		  
		  if (state == 4) //v6
		  {
			  for (int id_point=0; id_point< total_point; ++id_point) 
			  {
				  //int index = 76;
				  int index = 76+(8*id_point);
				  for (int ix=index; ix< index+2; ++ix) //76-84
				  {
					  if (read_buf[ix]>128)
					  {
						  signed char aaa;
						  //printf("ix now%d\n", ix);
						  aaa = read_buf[ix];
						  //printf("%d:\n", ix-index);
						  v6_point_before[ix-index] = aaa;
					  }
					  else
					  {
						  char aaa;
						  //printf("ix now%d\n", ix);
						  aaa = read_buf[ix];
						  //printf("%d:\n", ix-index);
						  v6_point_before[ix-index] = aaa;
					  }
				  }//printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8)));
				  for (int ix=index+2; ix< index+8; ix = ix+2) //76-84
				  {
					  if (((read_buf[ix] << 0) + (read_buf[ix+1] << 8))>32768)
					  {
						  signed short bbb;
						  //printf("ix now%d\n", ix);
						  bbb = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8));
						  v6_point_before[(ix-(index-2))/2] = bbb;
						  //printf("%d:\n", (ix-(index-2))/2);
					  }
					  else
					  {
						  short bbb;
						  //printf("ix now%d\n", ix);
						  bbb = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8));
						  v6_point_before[(ix-(index-2))/2] = bbb;
						  //printf("%d:\n", (ix-(index-2))/2);
					  }
				  }
				  size_t len_v6_point = sizeof(v6_point_before)/sizeof(v6_point_before[0]);
				  for(int num=0; num < len_v6_point; ++num)
				  {
					  v6_point_after[num] = v6_point_before[num]*V6_unit[num];
					  v6_2d_output[id_point][num] = v6_point_before[num]*V6_unit[num];
				  }
				  //printf("elevation:%f azimuth:%f doppler:%f range:%f snr:%f\n", v6_point_after[0], v6_point_after[1], v6_point_after[2], v6_point_after[3], v6_point_after[4]);				  
			      if (id_point == total_point - 1)
			      {
					  state = 5;		
				  }
			  }
		  }
		  if (state == 5)
		  {
			  FILE *fp = fopen(filename, "a");
			  //printf("state = %d\n", state);
			  int row = sizeof(v6_2d_output) / sizeof(v6_2d_output[0]);
			  int column = sizeof(v6_2d_output[0])/sizeof(v6_2d_output[0][0]);
			  int small_snr_count = 0;
			  float snr_tr = 2.0;
			  for(int num=0; num < row; ++num)//找snr小於2的
			  {
				  //v6_2d_output[num][4] ＝ v6_2d_output[num][4]-5.0;
				  if (v6_2d_output[num][4]<snr_tr)
				  {
					  printf("snr:%f\n", v6_2d_output[num][4]);
					  small_snr_count+=1;
				  }
			  }
			  int wo_snr = 	row - small_snr_count;
			  float v6_2d_output_wo_snr[wo_snr][5];
			  int count = 0;
			  for(int num=0; num < row; ++num)//把snr大於2的放陣列
			  {
				  if (v6_2d_output[num][4]>snr_tr)
				  {
					  v6_2d_output_wo_snr[count][0] = v6_2d_output[num][0];
					  v6_2d_output_wo_snr[count][1] = v6_2d_output[num][1];
					  v6_2d_output_wo_snr[count][2] = v6_2d_output[num][2];
					  v6_2d_output_wo_snr[count][3] = v6_2d_output[num][3];
					  v6_2d_output_wo_snr[count][4] = v6_2d_output[num][4];
					  count+=1;
				  }
			  } 
			  int row_wo_snr = sizeof(v6_2d_output_wo_snr) / sizeof(v6_2d_output_wo_snr[0]);
			  int column_wo_snr = sizeof(v6_2d_output_wo_snr[0])/sizeof(v6_2d_output_wo_snr[0][0]);

			  float pos1X[row_wo_snr][6];	 
			  float zOffSet = 1.0;
			  /*
			  printf("Number of rows: %d\n", row_wo_snr);
			  printf("Number of columns: %d\n", column_wo_snr);
			  */
			  for(int num=0; num < row_wo_snr; ++num)
			  {
				  float zt;
				  float xt;
				  float yt;
				  zt = v6_2d_output_wo_snr[num][3] * sin(v6_2d_output_wo_snr[num][0]) + zOffSet;
				  xt = v6_2d_output_wo_snr[num][3] * cos(v6_2d_output_wo_snr[num][0]) * sin(v6_2d_output_wo_snr[num][1]);
				  yt = v6_2d_output_wo_snr[num][3] * cos(v6_2d_output_wo_snr[num][0]) * cos(v6_2d_output_wo_snr[num][1]);
				  printf("elevation:%f azimuth:%f doppler:%f range:%f snr:%f\n", v6_2d_output[num][0], v6_2d_output[num][1], v6_2d_output[num][2], v6_2d_output[num][3], v6_2d_output[num][4]);				  
				  pos1X[num][0] = xt;
				  pos1X[num][1] = yt;
				  pos1X[num][2] = zt;
				  pos1X[num][3] = v6_2d_output_wo_snr[num][3];
				  pos1X[num][4] = v6_2d_output_wo_snr[num][2];
				  pos1X[num][5] = v6_2d_output_wo_snr[num][4];
				  //printf("x:%f \n", pos1X[num][0]);
				  //printf("y:%f \n", pos1X[num][1]);
				  //printf("z:%f \n", pos1X[num][2]);
				  //printf("range:%f \n", pos1X[num][3]);
				  //printf("Doppler:%f \n", pos1X[num][4]);
				  //printf("noise:%f \n", pos1X[num][5]);
				  //printf("x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5]);				  
				  //printf("\n");
			  }	
			  
			  output = dbscan_output(pos1X, row_wo_snr);

			  float index_point[row_wo_snr];
			  for(int num=0; num < row_wo_snr; ++num)
			  {
				  index_point[num] = output.vsos[num][3];
				  //printf("%f\n", output.vsos[num][3]);
			  }
			  
			  qsort(index_point, row_wo_snr, sizeof(float), cmpfunc);

			  int maxofindex = (int) index_point[row_wo_snr-1];
			  //printf("最大的數%d\n", maxofindex);
			  int countarray[maxofindex+1];
			  //陣列都設為0
			  for(int num=0; num < maxofindex+1; ++num)
			  {
				  countarray[num] = 0;
			  }
			  //找最多的dbscan標籤			  
			  for(int num=0; num < row_wo_snr; ++num)
			  {
				  //printf("dbscan的label:%f\n", index_point[num]);
				  for(int num1=0; num1 < maxofindex+1; ++num1)
				  {
					  //printf("num1 = %d\n", num1);
					  if ((int) index_point[num] == num1)
					  {
						  //printf("done\n");
						  countarray[num1]=countarray[num1]+1;
						  //printf("%d\n", countarray[num1]);
					  }
				  }
			  }
			  int max_index = 0;
			  for(int num=0; num < maxofindex+1; ++num)
			  {
				  if (countarray[num]>=countarray[max_index])
				  {
					  max_index = num;
				  }
			  }
			  //printf("最大的數:%d 總共有：%d\n", max_index, countarray[max_index]);
			  int max_index_num = countarray[max_index];
			  float sensorA[max_index_num][7];	
			  float x_array[max_index_num];
			  float y_array[max_index_num];
			  float z_array[max_index_num];
			  int index_sensorA = 0;
			  for(int num=0; num < row_wo_snr; ++num)
			  {
				  if ((int) output.vsos[num][3] == max_index)
				  {
				      //printf("x:%f y:%f z:%f range:%f Doppler:%f noise:%f label:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5], output.vsos[num][3]);
				      for(int num_1=0; num_1 < 6; ++num_1)
					  {
						  sensorA[index_sensorA][num_1] = pos1X[num][num_1];
					  }
					  x_array[index_sensorA] = pos1X[num][0];
					  y_array[index_sensorA] = pos1X[num][1];
					  z_array[index_sensorA] = pos1X[num][2];
					  sensorA[index_sensorA][6] = output.vsos[num][3];
					  index_sensorA+=1;		
				  }
			  }
			  
			  printf("index_sensorA %d\n", index_sensorA);
			  for(int num=0; num < max_index_num; ++num)
			  {
				   printf("x:%f y:%f z:%f range:%f Doppler:%f noise:%f label:%f\n", sensorA[num][0], sensorA[num][1], sensorA[num][2], sensorA[num][3], sensorA[num][4], sensorA[num][5], sensorA[num][6]);  
			  }
			  
			  qsort(x_array, max_index_num, sizeof(float), compare);
			  qsort(y_array, max_index_num, sizeof(float), compare);
			  qsort(z_array, max_index_num, sizeof(float), compare);
			  /*
			  for(int num=0; num < max_index_num; ++num)
			  {
				  printf("x_array = %f\n", x_array[num]);
			  }
			  */
			  float q1 = max_index_num*0.25;
			  float q3 = max_index_num*0.75;
			  printf("Q3:%f\n", q3);
			  printf("Q1:%f\n", q1);
			  int floor_q3, ceil_q3, floor_q1, ceil_q1;
			  floor_q3 = floor(q3);
			  ceil_q3 = ceil(q3);
			  printf("floor_q3:%d\n", floor_q3);
			  printf("ceil_q3:%d\n", ceil_q3);
			  floor_q1 = floor(q1);
			  ceil_q1 = ceil(q1);
			  printf("floor_q1:%d\n", floor_q1);
			  printf("ceil_q1:%d\n", ceil_q1);
			  float q1_x, q3_x, x_irq, x_min, x_max;
			  float q1_y, q3_y, y_irq, y_min, y_max;
			  float q1_z, q3_z, z_irq, z_min, z_max, z_mean, z_sum;
			  q1_x = (x_array[ceil_q1] + x_array[floor_q1])/2;
			  q3_x = (x_array[ceil_q3] + x_array[floor_q3])/2;
			  x_irq = q3_x - q1_x;
			  x_max = q3_x + 0.5*x_irq;
			  x_min = q1_x - 0.5*x_irq;
			  q1_y = (y_array[ceil_q1] + y_array[floor_q1])/2;
			  q3_y = (y_array[ceil_q3] + y_array[floor_q3])/2;
			  y_irq = q3_y - q1_y;
			  y_max = q3_y + 0.5*y_irq;
			  y_min = q1_y - 0.5*y_irq;
			  q1_z = (z_array[ceil_q1] + z_array[floor_q1])/2;
			  q3_z = (z_array[ceil_q3] + z_array[floor_q3])/2;
			  z_irq = q3_z - q1_z;
			  z_max = q3_z + 1.5*z_irq;
			  z_min = q1_z - 1.5*z_irq;
			  
			  printf("q1_x:%f\n", q1_x);
			  printf("q3_x:%f\n", q3_x);
			  printf("q1_y:%f\n", q1_y);
			  printf("q3_y:%f\n", q3_y);
			  printf("q1_z:%f\n", q1_z);
			  printf("q3_z:%f\n", q3_z);
			  printf("x_max = %f\n", x_max);
			  printf("x_min = %f\n", x_min);
			  printf("y_max = %f\n", y_max);
			  printf("y_min = %f\n", y_min);
			  printf("z_max = %f\n", z_max);
			  printf("z_min = %f\n", z_min);
			  
			  //find mean of array
			  z_sum = 0.0;
			  for(int num=0; num < max_index_num; ++num)
			  {
				  z_sum+=z_array[num]; 
			  }
			  z_mean = 0.0;
			  z_mean = z_sum/max_index_num;
			  float x_max_array[6], y_max_array[6], z_max_array[6], x_min_array[6], y_min_array[6], z_min_array[6], z_mean_array[6];
			  mean_count+=1;
			  //printf("mean_count = %d\n", mean_count);
			  float x_max_sum, x_min_sum, y_max_sum, y_min_sum, z_max_sum, z_min_sum, z_mean_sum;
			  if (mean_count>=5)
			  {   
				  //put number in last place
				  x_max_array[5] = x_max;
				  x_min_array[5] = x_min;
				  y_max_array[5] = y_max;
				  y_min_array[5] = y_min;
				  z_max_array[5] = z_max;
				  z_min_array[5] = z_min;
				  z_mean_array[5] = z_mean;
                  //shift
				  for(int num=0; num < 6; ++num)
				  {
					  x_max_array[num] = x_max_array[num+1];
					  x_min_array[num] = x_min_array[num+1];
					  y_max_array[num] = y_max_array[num+1];
					  y_min_array[num] = y_min_array[num+1];
					  z_max_array[num] = z_max_array[num+1];
					  z_min_array[num] = z_min_array[num+1];
					  z_mean_array[num] = z_mean_array[num+1];
				  }
				  x_max_sum = 0.0;
				  x_min_sum = 0.0;
				  y_max_sum = 0.0;
				  y_min_sum = 0.0;
				  z_max_sum = 0.0;
				  z_min_sum = 0.0;
				  z_mean_sum = 0.0;
				  for(int num=0; num < 5; ++num)
				  {
					  //printf("num = %d\n", num);
					  x_max_sum+=x_max_array[num]; 
					  x_min_sum+=x_min_array[num]; 
					  y_max_sum+=y_max_array[num]; 
					  y_min_sum+=y_min_array[num]; 
					  z_max_sum+=z_max_array[num]; 
					  z_min_sum+=z_min_array[num]; 
					  z_mean_sum+=z_mean_array[num]; 
				  }
				  x_max = x_max_sum/5;
				  x_min = x_min_sum/5;
				  y_max = y_max_sum/5;
				  y_min = y_min_sum/5;
				  z_max = z_max_sum/5;
				  z_min = z_min_sum/5;
				  z_mean = z_mean_sum/5;
				  /*
				  printf("x_max = %f\n", x_max);
				  printf("x_min = %f\n", x_min);
				  printf("y_max = %f\n", y_max);
				  printf("y_min = %f\n", y_min);
				  printf("z_max = %f\n", z_max);
				  printf("z_min = %f\n", z_min);
				  printf("z_mean = %f\n", z_mean);
				  */
				  //mean_count=0;
			  }
			  else
			  {
				  x_max_array[mean_count] = x_max;
				  x_min_array[mean_count] = x_min;
				  y_max_array[mean_count] = y_max;
				  y_min_array[mean_count] = y_min;
				  z_max_array[mean_count] = z_max;
				  z_min_array[mean_count] = z_min;
				  z_mean_array[mean_count] = z_mean;
			  }
			  float l1_dis;
			  l1_dis = sqrt(pow(x_min + ((x_max - x_min)/2), 2) + pow((y_min + (y_max - y_min)/2), 2));
			  float error_from_radar = 0.06;
			  z_mean = z_mean - error_from_radar;
			  float z_fall_lying[11];
			  if (fall_lying_count>10)
			  {
				  z_fall_lying[10] = z_mean;
				  for(int num=0; num < 10; ++num)
				  {
					  z_fall_lying[num] = z_fall_lying[num+1];
				  }
		      }
		      else
		      {
				  z_fall_lying[fall_lying_count] = z_mean;
			  }
			  int state_fall_lying = 0;
			  float tmp_threshold = 0.0;
			  float sum_z_fall_lying = 0.0;
			  float dif_z_thr = 0.0;
			  if (z_mean>0.7)
			  {
				  state_fall_lying = 0;
			  }
			  else
			  {
				  if (fall_lying_count>9)
				  {
					  for(int num=0; num < 5; ++num)
					  {
						  sum_z_fall_lying += z_fall_lying[num];
					  }
					  tmp_threshold = sum_z_fall_lying/5.0;
					  dif_z_thr = fabs(tmp_threshold-z_mean); 
					  if (dif_z_thr > 0.67)
					  {
						  state_fall_lying = 1;
					  }
				  }
			  }
			  /*
			  error_from_radar = randomForestModel.predict(np.array(l1_dis).reshape(-1, 1))
			  z_mean = z_mean - error_from_radar[0]
			  z_fall_lying.append(z_mean)
			  if len(z_fall_lying) > 10:
				  z_fall_lying.pop(0)
			  if z_mean > 0.8:  # z_men > 0.8 (reset to 0)
				  state_fall_lying = 0
			  else:
				  if len(z_fall_lying) > 9:
					  tmp_threshold = np.mean(z_fall_lying[:5])
					  dif_z_thr = np.abs(tmp_threshold - z_mean)
					  if dif_z_thr > 0.2:  # add: 放寬閾值
						  state_fall_lying = 1
			  */
			  float lenofx, lenofy, lenofz;
			  lenofz = z_max - z_min;
			  lenofx = x_max - x_min;
			  lenofy = y_max - y_min;
			  /*
			  printf("lenofz = %f\n", lenofz);
			  printf("lenofx = %f\n", lenofx);
			  printf("lenofy = %f\n", lenofy);
			  */
			  float stardand = 1.0;
			  int  state = 1;
			  
			  if (lenofz/lenofx >= 1.0 || lenofz/lenofy >= 1.0 || lenofx == 0.0 || lenofy == 0.0)
			  {
				  if (z_mean > (stardand / 10) * 7.3)
				  {
					  state = 1; //站
				  }
				  else if (lenofz/lenofx < 0.85 || lenofz/lenofy < 0.85)
				  {
					  if (state_fall_lying == 0)
					  {
						  state = 2; //臥
					  }
					  else
					  {
						  state = 3; //跌
					  }
				  }
				  else
				  {
					  state = 0; // 坐
				  }
			  }
			  else
			  {
				  if (state_fall_lying == 0)
				  {
					  state = 2; //臥
				  }
				  else
				  {
					  state = 3; //跌
				  }
			  }
			  //char state_str[4] = {"坐", "站", "臥", "跌"};
			  time(&rawtime);
			  info = localtime(&rawtime);
			  printf("姿態 = %d %s\n", state, asctime(info));
			  fprintf(fp, "%d, %s\n", state, asctime(info));
			  fclose(fp);
			  //*/
		  }

		}
  }
  


}
