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
#include "dbscan_animals.c"
#include <stdlib.h>
int frame_number = 0; 
int point_cnt_array[3] = {0};
//printf("%d%d%d", point_cnt_array[0], point_cnt_array[1], point_cnt_array[2]);
int row_temp = 0;
int frame_number_inf = 0; 
int animal_count = 0;
int temp_maxofindex = 0;
int temp_count_nan_normal = 0;
float temp_store_mean_xy [100][2];
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
			  printf("現在第%d frame!\n", frame_number_inf);
			  int row = sizeof(v6_2d_output) / sizeof(v6_2d_output[0]);
			  int column = sizeof(v6_2d_output[0])/sizeof(v6_2d_output[0][0]);
			  float v6_2d_output_bigdata[1000][5];	
			  
			  //printf("共有%d個點雲，現在frame：%d\n", row, frame_number);
			  /*
			  for(int num=0; num < row; ++num)//找snr小於2的
			  {
				  printf("EVERY FRAME: elevation:%f azimuth:%f doppler:%f range:%f snr:%f\n", v6_2d_output[num][0], v6_2d_output[num][1], v6_2d_output[num][2], v6_2d_output[num][3], v6_2d_output[num][4]);		
			  }
			  */
			  if (frame_number_inf>2)
			  {
				  int cnt_3 = point_cnt_array[1] + point_cnt_array[2] + row;
				  printf("總共%d個點雲\n", cnt_3);
				  float pos1a[cnt_3][5];
				  for(int num=0; num < point_cnt_array[1] + point_cnt_array[2]; ++num)
				  {
					  //printf("%d\n", num);
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  pos1a[num][num_1] = v6_2d_output_bigdata[num+point_cnt_array[0]][num_1];
					  }
				  }	
				  //printf("============");	
				  for(int num=0; num < row; ++num)
				  {
					  //printf("%d\n", num+point_cnt_array[1] + point_cnt_array[2]);
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  pos1a[num+point_cnt_array[1] + point_cnt_array[2]][num_1] = v6_2d_output[num][num_1];
					  }
				  }	
				  /*
				  for(int num=0; num < cnt_3; ++num)//找snr小於2的
				  {
					  printf("elevation:%f azimuth:%f doppler:%f range:%f snr:%f\n", pos1a[num][0], pos1a[num][1], pos1a[num][2], pos1a[num][3], pos1a[num][4]); //		
				  }	
				  */ 
				  int small_snr_count = 0;
				  float snr_tr = 8.0;
				  for(int num=0; num < cnt_3; ++num)//找snr小於2的
				  {
					  if (pos1a[num][4]<snr_tr)
					  {
						  printf("snr small detect!:%f\n", pos1a[num][4]);
						  small_snr_count+=1;
					  }
				  }
				  printf("small_snr_count%d\n", small_snr_count);
				  int wo_snr = 	cnt_3 - small_snr_count;
				  float pos1a_wo_snr[wo_snr][5];
				  int count = 0;
				  for(int num=0; num < cnt_3; ++num)//把snr大於2的放陣列
				  {
					  if (pos1a[num][4]>snr_tr)
					  {
						  pos1a_wo_snr[count][0] = pos1a[num][0];
						  pos1a_wo_snr[count][1] = pos1a[num][1];
						  pos1a_wo_snr[count][2] = pos1a[num][2];
						  pos1a_wo_snr[count][3] = pos1a[num][3];
						  pos1a_wo_snr[count][4] = pos1a[num][4];
						  count+=1;
					  }
				  } 
				  int row_wo_snr = sizeof(pos1a_wo_snr) / sizeof(pos1a_wo_snr[0]);
				  int column_wo_snr = sizeof(pos1a_wo_snr[0])/sizeof(pos1a_wo_snr[0][0]);
				  printf("row_wo_snr = %d\n", row_wo_snr);
				  float pos1X[row_wo_snr][6];	
				  int zero_nan_count = 0; 
				  for(int num=0; num < row_wo_snr; ++num)
				  {
					  float zt, yt, xt;
					  zt = 0.0;
					  xt = pos1a_wo_snr[num][3] * cos(pos1a_wo_snr[num][0]) * sin(pos1a_wo_snr[num][1]);
					  yt = pos1a_wo_snr[num][3] * cos(pos1a_wo_snr[num][0]) * cos(pos1a_wo_snr[num][1]);
					  pos1X[num][0] = xt;
					  pos1X[num][1] = yt;
					  pos1X[num][2] = zt;
					  pos1X[num][3] = pos1a_wo_snr[num][3];
					  pos1X[num][4] = pos1a_wo_snr[num][2];
					  pos1X[num][5] = pos1a_wo_snr[num][4];
					  //printf("x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5]);	
					  if ((pos1X[num][0]==0.0 && pos1X[num][1]==0.0 && pos1X[num][3]==0.0) || pos1X[num][0]== -0.0 || pos1X[num][1]== -0.0 || pos1X[num][3]== -0.0 || pos1X[num][4] < -10.0 || pos1X[num][4] > 10 || pos1X[num][0]+pos1X[num][1]>30.0 || pos1X[num][0]+pos1X[num][1]<-30.0)
					  {
						  zero_nan_count+=1;
					  }	  			  
				  }	
				  printf("zero_nan_count%d\n", zero_nan_count);
				  int wo_nan = row_wo_snr - zero_nan_count;
				  printf("wo_nan%d\n", wo_nan);
				  float pos1a_wo_nan[wo_nan][6];
				  int count_nan_normal = 0;
				  for(int num=0; num < row_wo_snr; ++num)//把snr大於2的放陣列
				  {
					  if ((pos1X[num][0]==0.0 && pos1X[num][1]==0.0 && pos1X[num][3]==0.0) || pos1X[num][0]== -0.0 || pos1X[num][1]== -0.0 || pos1X[num][3]== -0.0 || pos1X[num][4] < -10.0 || pos1X[num][4] > 10 || pos1X[num][0]+pos1X[num][1]>30.0 || pos1X[num][0]+pos1X[num][1]<-30.0)
					  {
						  //printf("pos1X[num][0] %f pos1X[num][1] %f pos1X[num][3] %f\n", pos1X[num][0], pos1X[num][1], pos1X[num][3]);
						  printf("DETECT!:x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5]);	
					  }
					  else
					  {
						  pos1a_wo_nan[count_nan_normal][0] = pos1X[num][0];
						  pos1a_wo_nan[count_nan_normal][1] = pos1X[num][1];
						  pos1a_wo_nan[count_nan_normal][2] = pos1X[num][2];
						  pos1a_wo_nan[count_nan_normal][3] = pos1X[num][3];
						  pos1a_wo_nan[count_nan_normal][4] = pos1X[num][4];
						  pos1a_wo_nan[count_nan_normal][5] = pos1X[num][5];
						  count_nan_normal+=1;						  
					  }
				  } 
				  
				  printf("pos1a_wo_nan\n");
				  for(int num=0; num < wo_nan; ++num)
				  {
					  printf("pos1a_wo_nan x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1a_wo_nan[num][0], pos1a_wo_nan[num][1], pos1a_wo_nan[num][2], pos1a_wo_nan[num][3], pos1a_wo_nan[num][4], pos1a_wo_nan[num][5]);	
				  }
				  
				  output = dbscan_output(pos1a_wo_nan, count_nan_normal);
				  
				  float index_point[wo_nan];
				  int dbscan_point[wo_nan];
				  //歸0陣列
				  for(int num=0; num < wo_nan; ++num)
				  {
					dbscan_point[num] = 0;  
				  }
				  //計算每個label有多少個點
				  for(int num=0; num < wo_nan; ++num)
				  {
					  index_point[num] = output.vsos[num][3];
					  for(int num_1=0; num_1 < wo_nan; ++num_1)
					  {
						  if (index_point[num] == num_1)
						  {
							  dbscan_point[num_1]+=1;
							  break;
						  }
					  }
				  }
				  
				  //把所有資料放在sensorA裡
				  float sensorA[wo_nan][7];	
				  for(int num=0; num < wo_nan; ++num)
				  {
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  sensorA[num][num_1] = pos1a_wo_nan[num][num_1];
					  }
					  sensorA[num][6] = index_point[num];
				  }
				  

				  animal_count+=1;
				  qsort(index_point, wo_nan, sizeof(float), cmpfunc);
				  int maxofindex = (int) index_point[wo_nan-1];
				  
				  printf("最大的數%d\n", maxofindex);
				  float store_mean_xy [maxofindex][2];
				  for(int num=0; num < maxofindex+1; ++num)
				  {
					  float mean_number_x = 0;
					  float mean_number_y = 0;
					  int label_count = 0;
					  for(int num_1=0; num_1 < wo_nan; ++num_1)
					  { 
						  if (num == sensorA[num_1][6])
						  {
							  printf("index=%d x=%f y=%f\n", num, sensorA[num_1][0],sensorA[num_1][1]);
							  //mean_number_x+=sensorA[num_1][0];
							  //mean_number_y+=sensorA[num_1][1]; 
							  label_count+=1;
						  }
					  }
					  float q_x[label_count];
					  float q_y[label_count];
					  int q_count = 0;
					  for(int num_1=0; num_1 < wo_nan; ++num_1)
					  { 
						  if (num == sensorA[num_1][6])
						  {
							  q_x[q_count]=sensorA[num_1][0];
							  q_y[q_count]=sensorA[num_1][1]; 
							  q_count+=1;
						  }
					  }
					  qsort(q_x, label_count, sizeof(float), compare);
					  qsort(q_y, label_count, sizeof(float), compare);
					  float q1 = label_count*0.25;
					  float q3 = label_count*0.75;
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

					  q1_x = (q_x[ceil_q1] + q_x[floor_q1])/2;
					  q3_x = (q_x[ceil_q3] + q_x[floor_q3])/2;
					  x_irq = q3_x - q1_x;
					  x_max = q3_x + 1.5*x_irq;
					  x_min = q1_x - 1.5*x_irq;
					  q1_y = (q_y[ceil_q1] + q_y[floor_q1])/2;
					  q3_y = (q_y[ceil_q3] + q_y[floor_q3])/2;
					  y_irq = q3_y - q1_y;
					  y_max = q3_y + 1.5*y_irq;
					  y_min = q1_y - 1.5*y_irq;
					  
					  printf("q1_x:%f\n", q1_x);
					  printf("q3_x:%f\n", q3_x);
					  printf("q1_y:%f\n", q1_y);
					  printf("q3_y:%f\n", q3_y);
					  printf("x_max = %f\n", x_max);
					  printf("x_min = %f\n", x_min);
					  printf("y_max = %f\n", y_max);
					  printf("y_min = %f\n", y_min);
					  int in_range_x, in_range_y;
					  for(int num_1=0; num_1 < wo_nan; ++num_1)
					  { 
						  if (num == sensorA[num_1][6])
						  {
							  if (sensorA[num_1][0]<x_max && sensorA[num_1][0]>x_min)
							  {
								  mean_number_x+=sensorA[num_1][0];
								  in_range_x+=1;
							  }
							  if (sensorA[num_1][1]<y_max && sensorA[num_1][1]>y_min)
							  {
								  mean_number_y+=sensorA[num_1][1];
								  in_range_y+=1;
							  }
							  
							  
						  }
					  }
					  //判斷nan
					  if (isnan(mean_number_x/in_range_x) == 1 || isnan(mean_number_y/in_range_y) == 1)
					  {
						  store_mean_xy[num][0] = 0.0; //
						  store_mean_xy[num][1] = 0.0;						  
					  }
					  else
					  {
						  store_mean_xy[num][0] = mean_number_x/in_range_x; //
						  store_mean_xy[num][1] = mean_number_y/in_range_y;						  
					  }
					  printf("==============\n");
				  }

				  
				  if (animal_count==1)
				  {
					  temp_maxofindex = (int) index_point[wo_nan-1];
					  temp_count_nan_normal = count_nan_normal;
					  
					  for(int num=0; num < maxofindex+1; ++num)
					  {
						  printf("index = %d mean x = %f mean_y = %f\n", num, store_mean_xy[num][0], store_mean_xy[num][1]);
					  }		
					  for(int num=0; num < temp_maxofindex+1; ++num)
					  {
						  temp_store_mean_xy[num][0] = store_mean_xy[num][0];
						  temp_store_mean_xy[num][1] = store_mean_xy[num][1];
					  }
				  }
				  else
				  {
					  printf("before frame = %d\n", temp_maxofindex);
					  printf("before num of point  = %d\n", temp_count_nan_normal);	
					  for(int num=0; num < temp_maxofindex+1; ++num)
					  {
						  printf("before index = %d mean x = %f mean_y = %f\n", num, temp_store_mean_xy[num][0], temp_store_mean_xy[num][1]);
					  }
					  printf("now frame = %d\n", maxofindex);
					  printf("now num of point  = %d\n", count_nan_normal);
					  for(int num=0; num < maxofindex+1; ++num)
					  {
						  printf("index = %d mean x = %f mean_y = %f\n", num, store_mean_xy[num][0], store_mean_xy[num][1]);
					  }						
					  if  (maxofindex == temp_maxofindex)
					  {
						  printf("cal dis:\n");
						  for(int num=0; num < maxofindex+1; ++num)
						  {
							  float temp_dis_x, temp_dis_y, dis;
							  temp_dis_x = pow((store_mean_xy[num][0] - temp_store_mean_xy[num][0]), 2);
							  temp_dis_y = pow((store_mean_xy[num][1] - temp_store_mean_xy[num][1]), 2);
							  dis = sqrt(temp_dis_x+temp_dis_y);
							  printf("index = %d, dis = %f\n", num, dis);
							  time(&rawtime);
							  info = localtime(&rawtime);							  
							  if (dis<=0.08)
							  {
								  printf("停止 %s\n", asctime(info));
								  fprintf(fp, "%d, 停止, %s\n", num, asctime(info));
							  }
							  else if(dis>=0.08 || dis<0.3)
							  {
								  printf("慢移 %s\n", asctime(info));
								  fprintf(fp, "%d, 慢移, %s\n", num, asctime(info));
							  } 
							  else
							  {
								  printf("快移 %s\n", asctime(info));
								  fprintf(fp, "%d, 快移, %s\n", num, asctime(info));
							  }
						  }
						  fclose(fp);
					  }
					  else
					  {
						  printf("else!\n");
						  for(int num=0; num < maxofindex+1; ++num)
						  {
							  printf("index = %d\n", num);
							  printf("x = %f y = %f\n", store_mean_xy[num][0], store_mean_xy[num][1]);
							  printf("慢移 %s\n", asctime(info));
							  fprintf(fp, "%d, 快移, %s\n", num, asctime(info));
						  }
						  fclose(fp);
						  
					  }
					  temp_maxofindex = maxofindex;		
					  temp_count_nan_normal = count_nan_normal;  
					  for(int num=0; num < 100; ++num)
					  {
						  temp_store_mean_xy[num][0] = 0.0;
						  temp_store_mean_xy[num][1] = 0.0;
					  }
					  for(int num=0; num < temp_maxofindex+1; ++num)
					  {
						  temp_store_mean_xy[num][0] = store_mean_xy[num][0];
						  temp_store_mean_xy[num][1] = store_mean_xy[num][1];
					  }
				  }

				  //---------------------------------------------
				  for(int num=0; num < 1000; ++num) // v6_2d_output_bigdata歸0
				  {
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  v6_2d_output_bigdata[num][num_1] = 0;
					  }
				  }		
				  for(int num=0; num < cnt_3; ++num) // pos1a放回v6_2d_output_bigdata
				  {
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  v6_2d_output_bigdata[num][num_1] = pos1a[num][num_1];
					  }
				  }
				  point_cnt_array[0] = point_cnt_array[1];
				  point_cnt_array[1] = point_cnt_array[2];
				  point_cnt_array[2] = row;

				  
			  }
			  else
			  {
				  point_cnt_array[frame_number] = row;
				  printf("row%d\n", row);
				  printf("row_temp%d\n", row_temp);
				  for(int num=0+row_temp; num < row+row_temp; ++num)
				  {
					  //printf("上一個frame的數量%d\n", row_temp);
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  v6_2d_output_bigdata[num][num_1] = v6_2d_output[num-row_temp][num_1];
					  }
					  //printf("第%d個frame, 共有%d個點雲\n",num, point_cnt_array[num]);		
				  }
				  row_temp = row+row_temp;
			  }
			  
			  for(int num=0; num < 3; ++num)//找snr小於2的
			  {
				  printf("第%d個frame, 共有%d個點雲\n",num, point_cnt_array[num]);		
			  }
			  if (frame_number == 3)
			  {
				  frame_number = 0;
			  }
			  frame_number += 1;
			  frame_number_inf+=1;
			  
		  }

		}
  }
}
