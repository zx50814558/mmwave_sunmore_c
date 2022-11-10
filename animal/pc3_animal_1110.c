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
int row_temp = 0;
int frame_number_inf = 0; 
int animal_count = 0;
int temp_maxofindex = 0;
int temp_count_nan_normal = 0;
int mode = 1; //顯示mode
float temp_store_mean_xy [100][2];
//int sort 的function
int compare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}
//float sort 的function
int cmpfunc (const void * a, const void * b)
{
	return (*(float*)a - *(float*)b);
}
//ieee 754 function decoder的strust 規範 ieee 754的格式
typedef union {
	float f;
	struct
	{
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} raw;
} myfloat;
//ieee 754 function 從二進位轉成int
unsigned int convertToInt(int* arr, int low, int high)
{
	unsigned f = 0, i;
	for (i = high; i >= low; i--) {
		f = f + arr[i] * pow(2, high - i);
	}
	return f;
}
// ieee 754 convert main function
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
// tlvTypeInfo的strust
struct tlvTypeInfo{
	int unitByte;
	int stateString;
	int sbyte;
	int dataByte;
	int retCnt;
	float nPoint;
	};
//玖邦python SDK寫死unitByte, sbyte, dataByte, stateString數值
struct tlvTypeInfo changeit(struct tlvTypeInfo, int state, int count);
struct tlvTypeInfo changeit(struct tlvTypeInfo s, int state, int count)
{
	s.unitByte = 20;
	s.sbyte = 8;
	s.dataByte = 0;
	s.stateString = 3;
	if (state = 2) //相當於python sdk v6
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
  //獲取時間
  time(&rawtime);
  info = localtime(&rawtime);
  //顯示時間
  printf("現在時間%s\n", asctime(info));
  char csv_name[10];
  printf("輸入檔名+.csv：");
  //輸入的檔名 變數=csv_name
  scanf("%s", csv_name);
  char *filename = csv_name;
  //宣告PORT號
  int serial_port = open("/dev/ttyTHS1", O_RDWR);
  struct termios tty;
  Struct output;
  //如果沒有PORT 跳錯誤
  if(tcgetattr(serial_port, &tty) != 0) {
      printf("Error %i from tcgetattr: %s\n", errno, strerror(errno));
      return 1;
  }
  //獲取時間
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
  //Baudrate
  cfsetispeed(&tty, B921600);
  cfsetospeed(&tty, B921600);
  //如果沒有PORT 跳錯誤
  if (tcsetattr(serial_port, TCSANOW, &tty) != 0) {
      printf("Error %i from tcsetattr: %s\n", errno, strerror(errno));
      return 1;
  }
  //放讀入的byte
  char read_buf [1024];
  int size;
  //限制一開始讀入的magicWord
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
		  if (mode == 0)
		  {
			printf("----------start-----------\n");
			printf("magic_words\n");
	      } 
		  //
		  //初始state設為0
		  int state = 0;
		  if (state == 0)
		  {
			  int same = 0;
			  //讀取8bytes
			  for (int ix=0; ix< 8; ++ix) 
			  {		
					//讀入的值要跟magicWord一樣才可以達到20hz
					if (read_buf[ix] == magicWord[ix])
					{   
						if (mode == 0)
						{
							printf("%d\n", read_buf[ix]);
						}
						same++;
						if (same==8)
						{   //讀滿8個magicWordstate轉成1
							state = 1;
						}
					}
			  }		  
		  }
		  
		  if (state ==1)
		  {	
			  if (mode == 0)
			  {
				  printf("----------header-----------\n");		
			  } 
			  //前8個已經被state=0讀過了，所以讀後40個bytes 1I=2H 1I=4Bytes
			  for (int ix=8; ix< 48; ++ix) 
			  {
				  if (ix<44)
				  {
				      //解self.hdr.version,self.hdr.totalPackLen,self.hdr.platform,self.hdr.frameNumber,self.hdr.subframeNumber,self.hdr.chirpMargin,self.hdr.frameMargin,self.hdr.trackProcessTime,self.hdr.uartSendTime,
					  if (ix%4==0)
					  {
						  if (mode == 0)
						  {
							  printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)) );
						  }
					  }	  
				  }
				  else
				  {
					  //解self.hdr.numTLVs,self.hdr.checksum
					  
					  if (ix%2==0)
					  {	  
						  if (mode == 0)
						  {
							  printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8)));
						  }
					  }  
					  if (ix==47)
					  {
						  state = 2;
						  //讀完40個byte後 state = 2
					  }
				  }
			  }
		  }
		  
		  if (state ==2)
		  {
			  int tlvLength;
			  if (mode ==0)
			  {
				  printf("----------TLV Header-----------\n");
			  }
			  for (int ix=48; ix< 56; ++ix) 
			  {
				  if (ix%4==0)
				  {
					  tlvLength = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)); 
					  if (mode == 0)
					  {
						  printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8) + (read_buf[ix+2] << 16) + (read_buf[ix+3] << 24)));
					  }
					  if (ix == 48)
					  {
						  //ttype 這裡是後來發現的這個值只會有6, 7, 8，對應sdk的state 
						  if (tlvLength != 6)
						  {
							  state = 0;
							  break;
						  }
					  }
					  if (ix == 52)
					  {
						  //tlvLength 超過10000會被重來
						  if (tlvLength > 10000)
						  {
							  state = 0;
							  break;
						  }
					  }
				  }
				  
				  if (ix==55)
				  {
					  //unitByteCount,lstate ,plen ,dataBytes,lenCount, numOfPoints = self.tlvTypeInfo(ttype,self.tlvLength,disp)丟進副程式tlvTypeInfo()
					  tlvTypeInfo_value1 = changeit(tlvTypeInfo_value, state, tlvLength);
					  // 會回傳以下值 供後續處理
					  if (mode == 0)
					  {
						  printf("unitByte:%d\n", tlvTypeInfo_value1.unitByte);
						  printf("stateString:%d\n", tlvTypeInfo_value1.stateString);
						  printf("sbyte:%d\n", tlvTypeInfo_value1.sbyte);
						  printf("dataByte:%d\n", tlvTypeInfo_value1.dataByte);
						  printf("retCnt:%d\n", tlvTypeInfo_value1.retCnt);
						  printf("nPoint:%f\n", tlvTypeInfo_value1.nPoint);						
					  }				  
					  state = tlvTypeInfo_value1.stateString;
				  }
				  
			  }
		  }
		  float V6_unit[5];
		  if (state == 3)
		  {
			  //self.u.elevationUnit,self.u.azimuthUnit,self.u.dopplerUnit,self.u.rangeUnit,self.u.snrUnit = struct.unpack('5f', sbuf) 解5個float
			  for (int ix2=0; ix2<5; ++ix2)  //5個float
			  {
				  for (int ix1=59+(ix2*4); ix1> 55+(ix2*4); --ix1) 
				  {
					  for (int ix=0; ix< 8; ++ix) 
					  {
						  ieee[k] = (read_buf[ix1] << ix & 0x80) >> 7; //把讀到的8個待處理的bytes位移存到陣列裡
						  k++;
					  }				  
				  }
				  //處理ieee 754
				  k = 0;
				  myfloat var;
				  unsigned f = convertToInt(ieee, 9, 31);
				  var.raw.mantissa = f;
				  f = convertToInt(ieee, 1, 8);      
				  var.raw.exponent = f;
				  var.raw.sign = ieee[0];
				  //V6_unit 存 self.u.elevationUnit,self.u.azimuthUnit,self.u.dopplerUnit,self.u.rangeUnit,self.u.snrUnit 
				  V6_unit[ix2] = var.f;
				  //printf("unit %f\n", var.f);				  
			  }	
			  state = 4;		  
		  }
		  int total_point;
		  total_point = (int) tlvTypeInfo_value1.retCnt/8.0;
		  int v6_point_before[5];
		  float v6_point_after[5];
		  float v6_2d_output[total_point][5];		  
		  //(e,a,d,r,s) = struct.unpack('2b3h', sbuf)   
		  if (state == 4) //v6
		  {
			
			  for (int id_point=0; id_point< total_point; ++id_point) 
			  {
				  //int index = 76;
				  int index = 76+(8*id_point);
				  for (int ix=index; ix< index+2; ++ix) //76-84 先解2b
				  {
					  if (read_buf[ix]>128)  //因為有正負號要先判斷
					  {
						  signed char aaa;
						  aaa = read_buf[ix];
						  v6_point_before[ix-index] = aaa;
					  }
					  else //因為有正負號要先判斷
					  {
						  char aaa;
						  aaa = read_buf[ix];
						  v6_point_before[ix-index] = aaa;
					  }
				  }//printf("%d\n", ((read_buf[ix] << 0) + (read_buf[ix+1] << 8)));
				  for (int ix=index+2; ix< index+8; ix = ix+2) //76-84 再解3h
				  {
					  if (((read_buf[ix] << 0) + (read_buf[ix+1] << 8))>32768) //因為有正負號要先判斷
					  {
						  signed short bbb;
						  bbb = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8));
						  v6_point_before[(ix-(index-2))/2] = bbb;
					  }
					  else //因為有正負號要先判斷
					  {
						  short bbb;
						  bbb = ((read_buf[ix] << 0) + (read_buf[ix+1] << 8));
						  v6_point_before[(ix-(index-2))/2] = bbb;
					  }
				  }
				  //v6_point_before放的就是e,a,d,r,s
				  /*
				  elv = e * self.u.elevationUnit
				  azi = a * self.u.azimuthUnit
		    	  dop = d * self.u.dopplerUnit
				  ran = r * self.u.rangeUnit
				  snr = s * self.u.snrUnit
				  */
				  size_t len_v6_point = sizeof(v6_point_before)/sizeof(v6_point_before[0]);
				  for(int num=0; num < len_v6_point; ++num)
				  {
					  v6_point_after[num] = v6_point_before[num]*V6_unit[num]; //相乘得到elv, azi, dop, ran, snr
					  v6_2d_output[id_point][num] = v6_point_before[num]*V6_unit[num]; //相乘得到elv, azi, dop, ran, snr
				  }
				  if (mode == 0)
				  {
					  printf("elevation:%f azimuth:%f doppler:%f range:%f snr:%f\n", v6_point_after[0], v6_point_after[1], v6_point_after[2], v6_point_after[3], v6_point_after[4]);
				  }		  
			      if (id_point == total_point - 1)
			      {
					  state = 5;		
				  }
			  }
		  }
		  if (state == 5)
		  {
			  if (mode == 0)
			  {
				  printf("現在第%d frame!\n", frame_number_inf);
			  }
			  
			  int row = sizeof(v6_2d_output) / sizeof(v6_2d_output[0]); //二維陣列的大小
			  int column = sizeof(v6_2d_output[0])/sizeof(v6_2d_output[0][0]); //二維陣列的大小
			  float v6_2d_output_bigdata[1000][5];	
			  
			  if (frame_number_inf>2) //前3frame會繼續累加點雲
			  {
				  int cnt_3 = point_cnt_array[1] + point_cnt_array[2] + row;
				  //printf("總共%d個點雲\n", cnt_3);
				  float pos1a[cnt_3][5];
				  //放前2 FRAME到最終陣列
				  for(int num=0; num < point_cnt_array[1] + point_cnt_array[2]; ++num)
				  {
					  //printf("%d\n", num);
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  pos1a[num][num_1] = v6_2d_output_bigdata[num+point_cnt_array[0]][num_1];
					  }
				  }	
				  //放第3FRAME到最終陣列
				  for(int num=0; num < row; ++num)
				  {
					  for(int num_1=0; num_1 < 5; ++num_1) // put number in array
					  {
						  pos1a[num+point_cnt_array[1] + point_cnt_array[2]][num_1] = v6_2d_output[num][num_1];
					  }
				  }	
				  int small_snr_count = 0;
				  float snr_tr = 8.0; //太小的SNR刪除閥值! 可調整
				  for(int num=0; num < cnt_3; ++num)//找snr小於8的
				  {
					  if (pos1a[num][4]<snr_tr)
					  {    
						  if (mode == 0)
						  {
							  printf("snr small detect!:%f\n", pos1a[num][4]); //小於8會被PRINT出來
						  }
						  small_snr_count+=1;
					  }
				  }
				  if (mode == 0)
				  {
					  printf("small_snr_count%d\n", small_snr_count); //陣列大小變了 因為小於2的要刪掉
				  }
				  int wo_snr = 	cnt_3 - small_snr_count;
				  float pos1a_wo_snr[wo_snr][5];
				  int count = 0;
				  for(int num=0; num < cnt_3; ++num)//把snr大於8的放陣列
				  {
					  if (pos1a[num][4]>snr_tr)
					  {
						  pos1a_wo_snr[count][0] = pos1a[num][0]; //再放一遍
						  pos1a_wo_snr[count][1] = pos1a[num][1]; //再放一遍
						  pos1a_wo_snr[count][2] = pos1a[num][2]; //再放一遍
						  pos1a_wo_snr[count][3] = pos1a[num][3]; //再放一遍
						  pos1a_wo_snr[count][4] = pos1a[num][4]; //再放一遍
						  count+=1;
					  }
				  } 
				  int row_wo_snr = sizeof(pos1a_wo_snr) / sizeof(pos1a_wo_snr[0]); //二維陣列的大小
				  int column_wo_snr = sizeof(pos1a_wo_snr[0])/sizeof(pos1a_wo_snr[0][0]); //二維陣列的大小
				  //printf("row_wo_snr = %d\n", row_wo_snr);
				  float pos1X[row_wo_snr][6];	
				  int zero_nan_count = 0; 
				  /*
				  printf("Number of rows: %d\n", row_wo_snr);
				  printf("Number of columns: %d\n", column_wo_snr);
				  Python Code
				  for i in range(len(pct)):
				  	  zt = pct[i][3] * np.sin(pct[i][0]) + zOffSet
					  xt = pct[i][3] * np.cos(pct[i][0]) * np.sin(pct[i][1])
					  yt = pct[i][3] * np.cos(pct[i][0]) * np.cos(pct[i][1])
					  pos1X[i] = (xt,yt,zt,pct[i][3],pct[i][2],pct[i][4]) # [x,y,z,range,Doppler,noise]
				  */
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
					  if (mode == 0)
					  {
						  printf("x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5]);	
					  }	
					  
					  //偵測0或是NAN或是INF
					  if ((pos1X[num][0]==0.0 && pos1X[num][1]==0.0 && pos1X[num][3]==0.0) || pos1X[num][0]== -0.0 || pos1X[num][1]== -0.0 || pos1X[num][3]== -0.0 || pos1X[num][4] < -10.0 || pos1X[num][4] > 10 || pos1X[num][0]+pos1X[num][1]>30.0 || pos1X[num][0]+pos1X[num][1]<-30.0)
					  {
						  zero_nan_count+=1;
					  }	  			  
				  }	
				  //偵測0或是NAN或是INF
				  //printf("zero_nan_count%d\n", zero_nan_count);
				  int wo_nan = row_wo_snr - zero_nan_count;
				  //printf("wo_nan%d\n", wo_nan);
				  float pos1a_wo_nan[wo_nan][6];
				  int count_nan_normal = 0;
				  for(int num=0; num < row_wo_snr; ++num)//偵測0或是NAN或是INF 不是的放新陣列
				  {
					  
					  if ((pos1X[num][0]==0.0 && pos1X[num][1]==0.0 && pos1X[num][3]==0.0) || pos1X[num][0]== -0.0 || pos1X[num][1]== -0.0 || pos1X[num][3]== -0.0 || pos1X[num][4] < -10.0 || pos1X[num][4] > 10 || pos1X[num][0]+pos1X[num][1]>30.0 || pos1X[num][0]+pos1X[num][1]<-30.0)
					  {
						  //printf("pos1X[num][0] %f pos1X[num][1] %f pos1X[num][3] %f\n", pos1X[num][0], pos1X[num][1], pos1X[num][3]);
					       if (mode == 0)
						   {
							  printf("DETECT!:x:%f y:%f z:%f range:%f Doppler:%f noise:%f\n", pos1X[num][0], pos1X[num][1], pos1X[num][2], pos1X[num][3], pos1X[num][4], pos1X[num][5]);
						   }
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
				  //送進dbscan
				  output = dbscan_output(pos1a_wo_nan, wo_nan);
				  if (mode ==0)
				  {
					  for(int num=0; num < wo_nan; ++num)
					  {
						  printf("num = %d", num);
						  printf("output dbscan x:%f y:%f z:%f index:%f\n", output.vsos[num][0], output.vsos[num][1], output.vsos[num][2], output.vsos[num][3]);	
					  }					
				  }

				  
				  float index_point[wo_nan];

				  //計算每個label有多少個點
				  for(int num=0; num < wo_nan; ++num)
				  {
					  index_point[num] = output.vsos[num][3];
				  }
				  //把所有資料放在sensorA裡
				  float sensorA[wo_nan][4];	
				  for(int num=0; num < wo_nan; ++num)
				  {
					  for(int num_1=0; num_1 < 3; ++num_1) // put number in array
					  {
						  sensorA[num][num_1] = output.vsos[num][num_1];
					  }
					  sensorA[num][3] = output.vsos[num][3];
				  }
				  if (mode == 0)
				  {
					for(int num=0; num < wo_nan; ++num)
					{
						printf("sensorA dbscan x:%f y:%f z:%f index:%f\n", sensorA[num][0], sensorA[num][1], sensorA[num][2], sensorA[num][3]);	
					}					
				  }

				  
				  animal_count+=1;
				  qsort(index_point, wo_nan, sizeof(float), cmpfunc);
				  int maxofindex = (int) index_point[wo_nan-1];
				  if (mode == 0)
				  {
					  printf("有%d個群\n\n", maxofindex);
				  }
				  float store_mean_xy [maxofindex][2];
				  for(int num=0; num < maxofindex+1; ++num)
				  {
					  float mean_number_x = 0;
					  float mean_number_y = 0;
					  int label_count = 0;
					  for(int num_2=0; num_2 < wo_nan; ++num_2)
					  { 
						  if (num == (int) sensorA[num_2][3])
						  {
							  if (sensorA[num_2][0]<15 && sensorA[num_2][1]<15 && sensorA[num_2][0] > -15 && sensorA[num_2][1] > -15)
							  {
								if (mode == 0)
								{
									printf("index=%d x=%f y=%f\n", num, sensorA[num_2][0],sensorA[num_2][1]);
								} 
								label_count+=1;
								mean_number_x += sensorA[num_2][0];
								mean_number_y += sensorA[num_2][1];
							  }
							  //找出相同LABEL有幾個

						  }
					  }
					  if (mode == 0)
					  {
						  printf("計算中心點中\n\n");
					  }
					  if (isnan(mean_number_x / label_count) == 1 || isnan(mean_number_y / label_count) == 1)
					  {
						  store_mean_xy[num][0] = 0.0; //如果異常就設0.0
						  store_mean_xy[num][1] = 0.0; //如果異常就設0.0				  
					  }
					  else
					  {
						  store_mean_xy[num][0] = mean_number_x / label_count; //找到中心點
						  store_mean_xy[num][1] = mean_number_y / label_count; //找到中心點
					  }
					  //printf("中心點 X = %f  ", store_mean_xy[num][0]);
					  //printf("中心點 Y = %f\n\n", store_mean_xy[num][1]);
					  /*
					  //創建當前LABRL數量的陣列
					  float q_x[label_count];
					  float q_y[label_count];
					  int q_count = 0;
					  for(int num_1=0; num_1 < wo_nan; ++num_1)
					  {   //只要是相同LABEL都放在 q_x Q_y
						  if (num == sensorA[num_1][6])
						  {
							  q_x[q_count]=sensorA[num_1][0];
							  q_y[q_count]=sensorA[num_1][1]; 
							  q_count+=1;
						  }
					  }
					  //對xy維作排列 算四分位數
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
					  x_max = q3_x + 1.5*x_irq; //算出x_max
					  x_min = q1_x - 1.5*x_irq; //算出x_min
					  q1_y = (q_y[ceil_q1] + q_y[floor_q1])/2;
					  q3_y = (q_y[ceil_q3] + q_y[floor_q3])/2;
					  y_irq = q3_y - q1_y;
					  y_max = q3_y + 1.5*y_irq; //算出y_max
					  y_min = q1_y - 1.5*y_irq; //算出y_min
					  
					  printf("q1_x:%f\n", q1_x);
					  printf("q3_x:%f\n", q3_x);
					  printf("q1_y:%f\n", q1_y);
					  printf("q3_y:%f\n", q3_y);
					  printf("x_max = %f\n", x_max);
					  printf("x_min = %f\n", x_min);
					  printf("y_max = %f\n", y_max);
					  printf("y_min = %f\n", y_min);
					  
					  int in_range_x, in_range_y;
					  //看有幾個在x_max, x_min, y_max, y_max裡面
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
						  store_mean_xy[num][0] = 0.0; //如果異常就設0.0
						  store_mean_xy[num][1] = 0.0; //如果異常就設0.0				  
					  }
					  else
					  {
						  store_mean_xy[num][0] = mean_number_x/in_range_x; //找到中心點
						  store_mean_xy[num][1] = mean_number_y/in_range_y;	//找到中心點	  
					  }

					  printf("========處理完中心點了！======\n");
					  */
				  }

				  
				  if (animal_count==1) //如果是第1 frame沒得比較
				  {
					  temp_maxofindex = (int) index_point[wo_nan-1];
					  temp_count_nan_normal = count_nan_normal;
					  /*
					  for(int num=0; num < maxofindex+1; ++num)
					  {
						  printf("index = %d mean x = %f mean_y = %f\n", num, store_mean_xy[num][0], store_mean_xy[num][1]);
					  }		
					  */
					  //temp_store_mean_xy就放上一frame的資料	
					  for(int num=0; num < temp_maxofindex+1; ++num)
					  {
						  temp_store_mean_xy[num][0] = store_mean_xy[num][0];
						  temp_store_mean_xy[num][1] = store_mean_xy[num][1];
					  }
				  }
				  else
				  {	  
					  if(mode == 0)
					  {
						printf("\n上一禎有%d個群\n", temp_maxofindex);
						printf("上一禎有%d個點\n", temp_count_nan_normal);	
						for(int num=0; num < temp_maxofindex+1; ++num)
						{   //顯示上一frame 資訊
							printf("上一禎 index = %d mean x = %f mean_y = %f\n", num, temp_store_mean_xy[num][0], temp_store_mean_xy[num][1]);
						}
						printf("\n現在禎有%d個群\n", maxofindex);
						printf("現在禎有%d個點\n", count_nan_normal);
							
						for(int num=0; num < maxofindex+1; ++num)
						{   //現在的資訊
							printf("現在禎 index = %d mean x = %f mean_y = %f\n", num, store_mean_xy[num][0], store_mean_xy[num][1]);
						}							
					  } 
					  if (mode == 1)
					  {
						printf("\e[1;1H");
						system("clear");							
					  }
	
					  if  (maxofindex == temp_maxofindex)
					  {
						  //printf("cal dis:\n"); //計算距離中
						  printf("==============================================\n");
						  printf("|                   總共%d個                  |\n", maxofindex+1);
						  printf("==============================================\n");
						  for(int num=0; num < maxofindex+1; ++num)
						  {
							  float temp_dis_x, temp_dis_y, dis;
							  temp_dis_x = pow((store_mean_xy[num][0] - temp_store_mean_xy[num][0]), 2);
							  temp_dis_y = pow((store_mean_xy[num][1] - temp_store_mean_xy[num][1]), 2);
							  dis = sqrt(temp_dis_x+temp_dis_y); //計算l1 dis
							  printf("|                  index = %d                 |\n", num+1);
							  printf("==============================================\n");
							  time(&rawtime);
							  info = localtime(&rawtime);		
							  //用l1 dis判斷狀態					  
							  if (dis<=0.04) //可調整
							  {
								  FILE *fp = fopen(filename, "a");
								  if (fp == NULL)
								  {
									  //printf("error");
									  return -1;
								  }
								  printf("|          停止 %s", asctime(info));
								  printf("==============================================\n");
								  fprintf(fp, "%d | %d, 停止, %s", num+1, maxofindex+1, asctime(info));
								  fclose(fp);
							  }
							  else if(dis>=0.04 || dis<0.15) //可調整
							  {
								  FILE *fp = fopen(filename, "a");
								  if (fp == NULL)
								  {
									  printf("error");
									  return -1;
								  }
								  printf("|          慢移 %s", asctime(info));
								  printf("==============================================\n");
								  fprintf(fp, "%d | %d, 慢移, %s", num+1, maxofindex+1, asctime(info));
								  fclose(fp);
							  } 
							  else
							  {
								  FILE *fp = fopen(filename, "a");
								  if (fp == NULL)
								  {
									  printf("error");
									  return -1;
								  }
								  printf("|          快移 %s", asctime(info));
								  printf("==============================================\n");
								  fprintf(fp, "%d | %d, 快移, %s", num+1, maxofindex+1, asctime(info));
								  fclose(fp);
							  }
						  }
						  FILE *fp = fopen(filename, "a");
						  if (fp == NULL)
						  {
						  	  printf("error");
							  return -1;
						  }
						  fprintf(fp, "end\n");
						  fclose(fp);
					  }
					  else
					  {
						  //printf("else!\n");
						  printf("==============================================\n");
						  printf("|                   總共%d個                  |\n", maxofindex+1);
						  printf("==============================================\n");
						  for(int num=0; num < maxofindex+1; ++num)
						  {
							  FILE *fp = fopen(filename, "a");
							  if (fp == NULL)
							  {
								  printf("error");
								  return -1;
							  }
							  printf("|                  index = %d                 |\n", num+1);
							  printf("==============================================\n");
							  //printf("x = %f y = %f\n", store_mean_xy[num][0], store_mean_xy[num][1]);
							  printf("|          慢移 %s", asctime(info));
							  printf("==============================================\n");
							  fprintf(fp, "%d | %d, 慢移, %s", num+1, maxofindex+1, asctime(info));
							  
							  fclose(fp);
						  }
						  FILE *fp = fopen(filename, "a");
						  if (fp == NULL)
						  {
						  	  printf("error");
							  return -1;
						  }
						  fprintf(fp, "end\n");
						  fclose(fp);
						  //printf("=====================結束=========================\n");
						  //printf("=====================結束=========================\n\n");
						  
						  
					  }
					  //全部處裡完之後 把store_mean_xy的點雲放到temp_store_mean_xy
					  temp_maxofindex = maxofindex;		
					  temp_count_nan_normal = count_nan_normal;  
					  for(int num=0; num < 100; ++num) //動物上限100
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
				  //point_cnt_array這個陣列是放點雲數量的
				  point_cnt_array[frame_number] = row;
				  //printf("row%d\n", row);
				  //printf("row_temp%d\n", row_temp);
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
			  /*
			  for(int num=0; num < 3; ++num)  //顯示累積3FRAME的點雲數量
			  {
				  printf("第%d個frame, 共有%d個點雲\n",num, point_cnt_array[num]);		
			  }
			  */
			  if (frame_number == 3) //重製frame_number
			  {
				  frame_number = 0;
			  }
			  frame_number += 1;
			  frame_number_inf+=1;
			  
		  }

		}
  }
}
