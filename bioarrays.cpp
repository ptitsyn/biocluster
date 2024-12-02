//Copyright(c) Andrey Ptitsyn, PBRC 2003
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>


#include "cluster.h"
#include "bioarrays.h"



static int grid_df1[19]={1,2,3,4,5,6,7,8,9,10,12,15,20,24,30,40,60,120,MAXINT};
static int grid_df2[34]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,60,120,MAXINT};
static float FDISTR005[19][34]={
161.4476, 199.5000, 215.7073, 224.5832, 230.1619, 233.9860, 236.7684, 238.8827, 240.5433, 241.8817, 243.9060, 245.9499, 248.0131, 249.0518, 250.0951, 251.1432, 252.1957, 253.2529, 254.3144,
18.5128, 19.0000, 19.1643, 19.2468, 19.2964, 19.3295, 19.3532, 19.3710, 19.3848, 19.3959, 19.4125, 19.4291, 19.4458, 19.4541, 19.4624, 19.4707, 19.4791, 19.4874, 19.4957,
10.1280, 9.5521, 9.2766, 9.1172, 9.0135, 8.9406, 8.8867, 8.8452, 8.8123, 8.7855, 8.7446, 8.7029, 8.6602, 8.6385, 8.6166, 8.5944, 8.5720, 8.5494, 8.5264,
7.7086, 6.9443, 6.5914, 6.3882, 6.2561, 6.1631, 6.0942, 6.0410, 5.9988, 5.9644, 5.9117, 5.8578, 5.8025, 5.7744, 5.7459, 5.7170, 5.6877, 5.6581, 5.6281,
6.6079, 5.7861, 5.4095, 5.1922, 5.0503, 4.9503, 4.8759, 4.8183, 4.7725, 4.7351, 4.6777, 4.6188, 4.5581, 4.5272, 4.4957, 4.4638, 4.4314, 4.3985, 4.3650,
5.9874, 5.1433, 4.7571, 4.5337, 4.3874, 4.2839, 4.2067, 4.1468, 4.0990, 4.0600, 3.9999, 3.9381, 3.8742, 3.8415, 3.8082, 3.7743, 3.7398, 3.7047, 3.6689,
5.5914, 4.7374, 4.3468, 4.1203, 3.9715, 3.8660, 3.7870, 3.7257, 3.6767, 3.6365, 3.5747, 3.5107, 3.4445, 3.4105, 3.3758, 3.3404, 3.3043, 3.2674, 3.2298,
5.3177, 4.4590, 4.0662, 3.8379, 3.6875, 3.5806, 3.5005, 3.4381, 3.3881, 3.3472, 3.2839, 3.2184, 3.1503, 3.1152, 3.0794, 3.0428, 3.0053, 2.9669, 2.9276,
5.1174, 4.2565, 3.8625, 3.6331, 3.4817, 3.3738, 3.2927, 3.2296, 3.1789, 3.1373, 3.0729, 3.0061, 2.9365, 2.9005, 2.8637, 2.8259, 2.7872, 2.7475, 2.7067,
4.9646, 4.1028, 3.7083, 3.4780, 3.3258, 3.2172, 3.1355, 3.0717, 3.0204, 2.9782, 2.9130, 2.8450, 2.7740, 2.7372, 2.6996, 2.6609, 2.6211, 2.5801, 2.5379,
4.8443, 3.9823, 3.5874, 3.3567, 3.2039, 3.0946, 3.0123, 2.9480, 2.8962, 2.8536, 2.7876, 2.7186, 2.6464, 2.6090, 2.5705, 2.5309, 2.4901, 2.4480, 2.4045,
4.7472, 3.8853, 3.4903, 3.2592, 3.1059, 2.9961, 2.9134, 2.8486, 2.7964, 2.7534, 2.6866, 2.6169, 2.5436, 2.5055, 2.4663, 2.4259, 2.3842, 2.3410, 2.2962,
4.6672, 3.8056, 3.4105, 3.1791, 3.0254, 2.9153, 2.8321, 2.7669, 2.7144, 2.6710, 2.6037, 2.5331, 2.4589, 2.4202, 2.3803, 2.3392, 2.2966, 2.2524, 2.2064,
4.6001, 3.7389, 3.3439, 3.1122, 2.9582, 2.8477, 2.7642, 2.6987, 2.6458, 2.6022, 2.5342, 2.4630, 2.3879, 2.3487, 2.3082, 2.2664, 2.2229, 2.1778, 2.1307,
4.5431, 3.6823, 3.2874, 3.0556, 2.9013, 2.7905, 2.7066, 2.6408, 2.5876, 2.5437, 2.4753, 2.4034, 2.3275, 2.2878, 2.2468, 2.2043, 2.1601, 2.1141, 2.0658,
4.4940, 3.6337, 3.2389, 3.0069, 2.8524, 2.7413, 2.6572, 2.5911, 2.5377, 2.4935, 2.4247, 2.3522, 2.2756, 2.2354, 2.1938, 2.1507, 2.1058, 2.0589, 2.0096,
4.4513, 3.5915, 3.1968, 2.9647, 2.8100, 2.6987, 2.6143, 2.5480, 2.4943, 2.4499, 2.3807, 2.3077, 2.2304, 2.1898, 2.1477, 2.1040, 2.0584, 2.0107, 1.9604,
4.4139, 3.5546, 3.1599, 2.9277, 2.7729, 2.6613, 2.5767, 2.5102, 2.4563, 2.4117, 2.3421, 2.2686, 2.1906, 2.1497, 2.1071, 2.0629, 2.0166, 1.9681, 1.9168,
4.3807, 3.5219, 3.1274, 2.8951, 2.7401, 2.6283, 2.5435, 2.4768, 2.4227, 2.3779, 2.3080, 2.2341, 2.1555, 2.1141, 2.0712, 2.0264, 1.9795, 1.9302, 1.8780,
4.3512, 3.4928, 3.0984, 2.8661, 2.7109, 2.5990, 2.5140, 2.4471, 2.3928, 2.3479, 2.2776, 2.2033, 2.1242, 2.0825, 2.0391, 1.9938, 1.9464, 1.8963, 1.8432,
4.3248, 3.4668, 3.0725, 2.8401, 2.6848, 2.5727, 2.4876, 2.4205, 2.3660, 2.3210, 2.2504, 2.1757, 2.0960, 2.0540, 2.0102, 1.9645, 1.9165, 1.8657, 1.8117,
4.3009, 3.4434, 3.0491, 2.8167, 2.6613, 2.5491, 2.4638, 2.3965, 2.3419, 2.2967, 2.2258, 2.1508, 2.0707, 2.0283, 1.9842, 1.9380, 1.8894, 1.8380, 1.7831,
4.2793, 3.4221, 3.0280, 2.7955, 2.6400, 2.5277, 2.4422, 2.3748, 2.3201, 2.2747, 2.2036, 2.1282, 2.0476, 2.0050, 1.9605, 1.9139, 1.8648, 1.8128, 1.7570,
4.2597, 3.4028, 3.0088, 2.7763, 2.6207, 2.5082, 2.4226, 2.3551, 2.3002, 2.2547, 2.1834, 2.1077, 2.0267, 1.9838, 1.9390, 1.8920, 1.8424, 1.7896, 1.7330,
4.2417, 3.3852, 2.9912, 2.7587, 2.6030, 2.4904, 2.4047, 2.3371, 2.2821, 2.2365, 2.1649, 2.0889, 2.0075, 1.9643, 1.9192, 1.8718, 1.8217, 1.7684, 1.7110,
4.2252, 3.3690, 2.9752, 2.7426, 2.5868, 2.4741, 2.3883, 2.3205, 2.2655, 2.2197, 2.1479, 2.0716, 1.9898, 1.9464, 1.9010, 1.8533, 1.8027, 1.7488, 1.6906,
4.2100, 3.3541, 2.9604, 2.7278, 2.5719, 2.4591, 2.3732, 2.3053, 2.2501, 2.2043, 2.1323, 2.0558, 1.9736, 1.9299, 1.8842, 1.8361, 1.7851, 1.7306, 1.6717,
4.1960, 3.3404, 2.9467, 2.7141, 2.5581, 2.4453, 2.3593, 2.2913, 2.2360, 2.1900, 2.1179, 2.0411, 1.9586, 1.9147, 1.8687, 1.8203, 1.7689, 1.7138, 1.6541,
4.1830, 3.3277, 2.9340, 2.7014, 2.5454, 2.4324, 2.3463, 2.2783, 2.2229, 2.1768, 2.1045, 2.0275, 1.9446, 1.9005, 1.8543, 1.8055, 1.7537, 1.6981, 1.6376,
4.1709, 3.3158, 2.9223, 2.6896, 2.5336, 2.4205, 2.3343, 2.2662, 2.2107, 2.1646, 2.0921, 2.0148, 1.9317, 1.8874, 1.8409, 1.7918, 1.7396, 1.6835, 1.6223,
4.0847, 3.2317, 2.8387, 2.6060, 2.4495, 2.3359, 2.2490, 2.1802, 2.1240, 2.0772, 2.0035, 1.9245, 1.8389, 1.7929, 1.7444, 1.6928, 1.6373, 1.5766, 1.5089,
4.0012, 3.1504, 2.7581, 2.5252, 2.3683, 2.2541, 2.1665, 2.0970, 2.0401, 1.9926, 1.9174, 1.8364, 1.7480, 1.7001, 1.6491, 1.5943, 1.5343, 1.4673, 1.3893,
3.9201, 3.0718, 2.6802, 2.4472, 2.2899, 2.1750, 2.0868, 2.0164, 1.9588, 1.9105, 1.8337, 1.7505, 1.6587, 1.6084, 1.5543, 1.4952, 1.4290, 1.3519, 1.2539,
3.8415, 2.9957, 2.6049, 2.3719, 2.2141, 2.0986, 2.0096, 1.9384, 1.8799, 1.8307, 1.7522, 1.6664, 1.5705, 1.5173, 1.4591, 1.3940, 1.3180, 1.2214, 1.0000
};
////////////////////////////////////////////////////////////////////////////


int float_sorter( const void *x, const void *y)
{
	return((int)(*(float *)x-*(float *)y));
}

int float_pointer_sorter( void **x, void **y)
{
	return((int)(**(float **)x-**(float **)y));
}

int int_sorter( int *x, int *y) {
	return((int)(*x-*y));
}

char *FGetTabS(FILE *f){
int i;
char c;
char s[1000];
char *str;
	i=0;
	while(1){
		fscanf(f,"%c",&c);
		s[i]=0;
      if(c=='\t') break;
      s[i++]=c;
   }
   str=(char *)malloc(strlen(s)+1);
   strcpy(str,s);
   return(str);
};

char *GetTabS(char *s){
char buf[1000];
char *ss;
char *str;
	strcpy(buf,s);
   ss=strchr(buf,'\t');
	if(ss){
   	*ss=0;
   }else{
	   ss=strchr(buf,'\n');
      if(ss){
      	*ss=0;
		}else{
	   	return(NULL);
      }
   }
   str=(char *)malloc(strlen(buf)+1);
   strcpy(str,buf);
   return(str);
};

float F_table( int x, int y){
int i,j;
	i=j=0;
   while(x>=grid_df1[i]) i++;
   while(y>=grid_df2[j]) j++;
   return(FDISTR005[i][j]);
};


//////////////////////////////////////////////////////////////////////////
ScanArrayHeader::ScanArrayHeader(){
};

ScanArrayHeader::~ScanArrayHeader(){
	if(UserName) delete UserName;
	if(Computer) delete Computer;
	if(Date) delete Date;
	if(Experiment) delete Experiment;
	if(ExperimentPath) delete ExperimentPath;
	if(Protocol) delete Protocol;
};

ScanArrayHeader::Find(FILE *f){
char s[1000];
   rewind(f);
  	if(!fgets(s,1000,f)){
		return(0);
   }
   if(strstr(s,"User Name")) {
		rewind(f);
   	return(1);
	}else{
	   return(0);
   }
};

ScanArrayHeader::Read(FILE *f){
	return(0);
};

ScanArrayHeader::Write(FILE *f){
	return(0);
};

/////////////////////////////////////////////////////////////////////////

ProtocolInfo::ProtocolInfo(){
};

ProtocolInfo::~ProtocolInfo(){
	if(Units) delete Units;
   if(InterstitialComment) delete InterstitialComment;
	if(QuantificationMethod) delete QuantificationMethod;
	if(QualityConfidenceCalculation) delete QualityConfidenceCalculation;
};

ProtocolInfo::Find(FILE *f){
char s[1000];
int flag;
	flag=0;
	while(flag<2){
   	if(!fgets(s,1000,f)){
      	rewind(f);
         flag++;
      }
      if(strstr(s,"Begin Protocol Info")) return(flag);
	}
   return(flag);
};

char *ffgets(char *s, int l, FILE *f){
	while(1){
   	if(!fgets(s,l,f)) return(NULL);
   	if(isalpha(s[0])) return(s);
   }
}

int ProtocolInfo::Read(FILE *f){
char s[1000];
	while(1){
		if(!fgets(s,1000,f)) return(0);
      if(strstr(s,"End Protocol")) break;
  		if(strstr(s,"Units")) {
      	Units=GetTabS(strchr(s,'\t')+1);
         continue;
      }
  		if(strstr(s,"Array Columns Spacing")){
      	sscanf(strchr(s,'\t'),"%f", &ArrayColumnSpacing);
         continue;
      }
  		if(strstr(s,"Array Row Spacing")){
      	sscanf(strchr(s,'\t'),"%f", &ArrayRowSpacing);
         continue;
      }
  		if(strstr(s,"Spot Rows Spacing")){
      	sscanf(strchr(s,'\t'),"%d", &SpotRowSpacing);
         continue;
      }
	  	if(strstr(s,"Spot Columns Spacing")){
      	sscanf(strchr(s,'\t'),"%d", &SpotColumnSpacing);
         continue;
      }
  		if(strstr(s,"Array Rows")) {
      	sscanf(strchr(s,'\t'),"%d", &ArrayRows);
         continue;
      }
	  	if(strstr(s,"Array Columns")) {
      	sscanf(strchr(s,'\t'),"%d", &ArrayColumns);
         continue;
      }
	  	if(strstr(s,"Rows")){
      	sscanf(strchr(s,'\t'),"%d", &Rows);
         continue;
      }
  		if(strstr(s,"Columns")){
      	sscanf(strchr(s,'\t'),"%d", &Columns);
         continue;
      }
  		if(strstr(s,"Spot Diameter")){
      	sscanf(strchr(s,'\t'),"%d", &SpotDiameter);
         continue;
      }
	  	if(strstr(s,"Interstitial")){
      	sscanf(strchr(s,'\t'),"%d", &Interstitial);
	  		InterstitialComment=GetTabS(strchr(s,'\"')+1);
		   *strchr(InterstitialComment,'\"')=0;
         continue;
      }
	  	if(strstr(s,"Spots Per")) {
   		sscanf(strchr(s,'\t'),"%d", &SpotsPerArray);
      	continue;
	   }
  		if(strstr(s,"Total Spots")){
   		sscanf(strchr(s,'\t'),"%d", &TotalSpots);
	      continue;
   	}
	  	if(strstr(s,"Crosstalk")) {
	   	if(!strstr(s,"not")) CrosstalkCorrected=1;
	      continue;
	   }
  		if(strstr(s,"Background")) {
	   	if (!strstr(s,"not")) BackgroundSubtracted=1;
	      continue;
      }
	 	if(strstr(s,"Quantification")){
   		QuantificationMethod=GetTabS(strchr(s,'\t')+1);
      	continue;
	   }
	  	if(strstr(s,"Quality Confidence")) {
   		QualityConfidenceCalculation=GetTabS(strchr(s,'\t')+1);
      	continue;
	   }
   }
  	return(1);
};

int ProtocolInfo::Write(FILE *f){
	return(0);
}

////////////////////////////////////////////////////////////////////////
ToleranceWeight::ToleranceWeight(){
};

ToleranceWeight::~ToleranceWeight(){
};

int ToleranceWeight::Find(FILE *){
	return(0);
};

int ToleranceWeight::Read(FILE *){
	return(0);
};

int ToleranceWeight::Write(FILE *){
	return(0);
};

/////////////////////////////////////////////////////////////////////

ImageInfo::ImageInfo(){
};

ImageInfo::~ImageInfo(){
};

ImageInfo::Find(FILE *f){
};

ImageInfo::Read(FILE *f){
};

ImageInfo::Write(FILE *f){
};

//////////////////////////////////////////////////////////////////////

ScanArraySpot::ScanArraySpot(){
	IgnoreFilter=1;
};

ScanArraySpot::~ScanArraySpot(){
	if(Name) delete(Name);
};


ScanArraySpot::Read(FILE *f, int num){
char ss[4096];
char *s,*sp;
int i,j;
	s=ss;
	fgets(s,1024,f);
	i=0;
	i+=sscanf(s,"%d",&Number);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%d\t",&ArrayRow);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%d\t",&ArrayColumn);
   s=strchr(s,'\t')+1;
	i+=sscanf(s+1,"%d\t",&Row);
   s=strchr(s,'\t')+1;
	i+=sscanf(s+1,"%d\t",&Column);
   s=strchr(s,'\t')+1;
   Name=(char *)calloc(strchr(s,'\t')-s+1,sizeof(char));
	memcpy(Name,s,strchr(s,'\t')-s);
	i++;
   s=strchr(s,'\t')+1;
	//Name=FGetTabS(f);
	i+=sscanf(s,"%d\t",&XLocation);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%d\t",&YLocation);
   s=strchr(s,'\t')+1;
//   if(i!=7) return(IgnoreFilter=1);

	i=0;
	i+=sscanf(s,"%f\t",&ch1Intensity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1Background);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1IntensityStdDev);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1BackgroundStdDev);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1Diameter);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%d\t",&ch1Area);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1Footprint);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1Circularity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1SpotUniformity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1BkgUniformity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1SignalNoiseRatio);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch1Confidence);
   s=strchr(s,'\t')+1;
   if(i!=12) return(IgnoreFilter=1);

	i=0;
	i+=sscanf(s,"%f\t",&ch2Intensity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2Background);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2IntensityStdDev);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2BackgroundStdDev);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2Diameter);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%d\t",&ch2Area);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2Footprint);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2Circularity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2SpotUniformity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2BkgUniformity);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2SignalNoiseRatio);
   s=strchr(s,'\t')+1;
	i+=sscanf(s,"%f\t",&ch2Confidence);
   s=strchr(s,'\t')+1;
   if(i!=12) return(IgnoreFilter=1);

	sp=s;
   while(sp=strchr(sp,'\t')){;
		s=sp;
   	sp++;
   }
	sscanf(s,"%d",&IgnoreFilter);
   //s=strchr(s,'\t')+1;
   IgnoreFilter=(short int)((IgnoreFilter+1)%2);
   return(IgnoreFilter);

/*
	if(num==33){
		sscanf(s,"%d\t",&IgnoreFilter);
	   s=strchr(s,'\t')+1;
	   IgnoreFilter=(short int)((IgnoreFilter+1)%2);
      return(IgnoreFilter);
   }else{
		i=0;
		i+=sscanf(s,"%d\t",&ch1DiameterFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1AreaFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1FootprintFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1CircularityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1SpotUniformityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1BkgUniformityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1SignalNoiseRatioFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch1ReplicateFilter);
	   s=strchr(s,'\t')+1;

		i+=sscanf(s,"%d\t",&ch2DiameterFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2AreaFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2FootprintFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2CircularityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2SpotUniformityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2BkgUniformityFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2SignalNoiseRatioFilter);
	   s=strchr(s,'\t')+1;
		i+=sscanf(s,"%d\t",&ch2ReplicateFilter);
	   s=strchr(s,'\t')+1;

		i+=sscanf(s,"%d\t",&IgnoreFilter);
   }

   if(i<17) return(IgnoreFilter=1);
   ch1DiameterFilter=(short int)((ch1DiameterFilter+1)%2);
   ch1AreaFilter=(short int)((ch1AreaFilter+1)%2);
   ch1FootprintFilter=(short int)((ch1FootprintFilter+1)%2);
   ch1CircularityFilter=(short int)((ch1CircularityFilter+1)%2);
   ch1SpotUniformityFilter=(short int)((ch1SpotUniformityFilter+1)%2);
   ch1BkgUniformityFilter=(short int)((ch1BkgUniformityFilter+1)%2);
   ch1SignalNoiseRatioFilter=(short int)((ch1SignalNoiseRatioFilter+1)%2);
   ch1ReplicateFilter=(short int)((ch1ReplicateFilter+1)%2);
   ch2DiameterFilter=(short int)((ch2DiameterFilter+1)%2);
   ch2AreaFilter=(short int)((ch2AreaFilter+1)%2);
   ch2FootprintFilter=(short int)((ch2FootprintFilter+1)%2);
   ch2CircularityFilter=(short int)((ch2CircularityFilter+1)%2);
   ch2SpotUniformityFilter=(short int)((ch2SpotUniformityFilter+1)%2);
   ch2BkgUniformityFilter=(short int)((ch2BkgUniformityFilter+1)%2);
   ch2SignalNoiseRatioFilter=(short int)((ch2SignalNoiseRatioFilter+1)%2);
   ch2ReplicateFilter=(short int)((ch2ReplicateFilter+1)%2);
   IgnoreFilter=(short int)((IgnoreFilter+1)%2);
   */
};

ScanArraySpot::Write(FILE *f){
	fprintf(f,"%d\t",Number);
	fprintf(f,"%d\t",ArrayRow);
	fprintf(f,"%d\t",ArrayColumn);
	fprintf(f,"%d\t",Row);
	fprintf(f,"%d\t",Column);
	fprintf(f,"%s\t",Name);
	fprintf(f,"%d\t",XLocation);
	fprintf(f,"%d\t",YLocation);

	fprintf(f,"%f\t",ch1Intensity);
	fprintf(f,"%f\t",ch1Background);
	fprintf(f,"%f\t",ch1IntensityStdDev);
	fprintf(f,"%f\t",ch1BackgroundStdDev);
	fprintf(f,"%f\t",ch1Diameter);
	fprintf(f,"%d\t",ch1Area);
	fprintf(f,"%f\t",ch1Footprint);
	fprintf(f,"%f\t",ch1Circularity);
	fprintf(f,"%f\t",ch1SpotUniformity);
	fprintf(f,"%f\t",ch1BkgUniformity);
	fprintf(f,"%f\t",ch1SignalNoiseRatio);
	fprintf(f,"%f\t",ch1Confidence);

	fprintf(f,"%f\t",ch2Intensity);
	fprintf(f,"%f\t",ch2Background);
	fprintf(f,"%f\t",ch2IntensityStdDev);
	fprintf(f,"%f\t",ch2BackgroundStdDev);
	fprintf(f,"%f\t",ch2Diameter);
	fprintf(f,"%d\t",ch2Area);
	fprintf(f,"%f\t",ch2Footprint);
	fprintf(f,"%f\t",ch2Circularity);
	fprintf(f,"%f\t",ch2SpotUniformity);
	fprintf(f,"%f\t",ch2BkgUniformity);
	fprintf(f,"%f\t",ch2SignalNoiseRatio);
	fprintf(f,"%f\t",ch2Confidence);

	fprintf(f,"%d\t",ch1DiameterFilter);
	fprintf(f,"%d\t",ch1AreaFilter);
	fprintf(f,"%d\t",ch1FootprintFilter);
	fprintf(f,"%d\t",ch1CircularityFilter);
	fprintf(f,"%d\t",ch1SpotUniformityFilter);
	fprintf(f,"%d\t",ch1BkgUniformityFilter);
	fprintf(f,"%d\t",ch1SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch1ReplicateFilter);

	fprintf(f,"%d\t",ch2DiameterFilter);
	fprintf(f,"%d\t",ch2AreaFilter);
	fprintf(f,"%d\t",ch2FootprintFilter);
	fprintf(f,"%d\t",ch2CircularityFilter);
	fprintf(f,"%d\t",ch2SpotUniformityFilter);
	fprintf(f,"%d\t",ch2BkgUniformityFilter);
	fprintf(f,"%d\t",ch2SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch2ReplicateFilter);

	fprintf(f,"%d\n",IgnoreFilter);

	return(IgnoreFilter);
};

ScanArraySpot::Write(FILE *f, float extravalue){
	fprintf(f,"%d\t",Number);
	fprintf(f,"%d\t",ArrayRow);
	fprintf(f,"%d\t",ArrayColumn);
	fprintf(f,"%d\t",Row);
	fprintf(f,"%d\t",Column);
	fprintf(f,"%s\t",Name);
	fprintf(f,"%d\t",XLocation);
	fprintf(f,"%d\t",YLocation);

   fprintf(f,"%f\t", extravalue);

	fprintf(f,"%f\t",ch1Intensity);
	fprintf(f,"%f\t",ch1Background);
	fprintf(f,"%f\t",ch1IntensityStdDev);
	fprintf(f,"%f\t",ch1BackgroundStdDev);
	fprintf(f,"%f\t",ch1Diameter);
	fprintf(f,"%d\t",ch1Area);
	fprintf(f,"%f\t",ch1Footprint);
	fprintf(f,"%f\t",ch1Circularity);
	fprintf(f,"%f\t",ch1SpotUniformity);
	fprintf(f,"%f\t",ch1BkgUniformity);
	fprintf(f,"%f\t",ch1SignalNoiseRatio);
	fprintf(f,"%f\t",ch1Confidence);

	fprintf(f,"%f\t",ch2Intensity);
	fprintf(f,"%f\t",ch2Background);
	fprintf(f,"%f\t",ch2IntensityStdDev);
	fprintf(f,"%f\t",ch2BackgroundStdDev);
	fprintf(f,"%f\t",ch2Diameter);
	fprintf(f,"%d\t",ch2Area);
	fprintf(f,"%f\t",ch2Footprint);
	fprintf(f,"%f\t",ch2Circularity);
	fprintf(f,"%f\t",ch2SpotUniformity);
	fprintf(f,"%f\t",ch2BkgUniformity);
	fprintf(f,"%f\t",ch2SignalNoiseRatio);
	fprintf(f,"%f\t",ch2Confidence);

	fprintf(f,"%d\t",ch1DiameterFilter);
	fprintf(f,"%d\t",ch1AreaFilter);
	fprintf(f,"%d\t",ch1FootprintFilter);
	fprintf(f,"%d\t",ch1CircularityFilter);
	fprintf(f,"%d\t",ch1SpotUniformityFilter);
	fprintf(f,"%d\t",ch1BkgUniformityFilter);
	fprintf(f,"%d\t",ch1SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch1ReplicateFilter);

	fprintf(f,"%d\t",ch2DiameterFilter);
	fprintf(f,"%d\t",ch2AreaFilter);
	fprintf(f,"%d\t",ch2FootprintFilter);
	fprintf(f,"%d\t",ch2CircularityFilter);
	fprintf(f,"%d\t",ch2SpotUniformityFilter);
	fprintf(f,"%d\t",ch2BkgUniformityFilter);
	fprintf(f,"%d\t",ch2SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch2ReplicateFilter);

	fprintf(f,"%d\n",IgnoreFilter);

	return(IgnoreFilter);
};

ScanArraySpot::Write(FILE *f, float extravalue1, float extravalue2, float extravalue3){
	fprintf(f,"%d\t",Number);
	fprintf(f,"%d\t",ArrayRow);
	fprintf(f,"%d\t",ArrayColumn);
	fprintf(f,"%d\t",Row);
	fprintf(f,"%d\t",Column);
	fprintf(f,"%s\t",Name);
	fprintf(f,"%d\t",XLocation);
	fprintf(f,"%d\t",YLocation);

   fprintf(f,"%f\t", extravalue1);
   fprintf(f,"%f\t", extravalue2);
   fprintf(f,"%f\t", extravalue3);

	fprintf(f,"%f\t",ch1Intensity);
	fprintf(f,"%f\t",ch1Background);
	fprintf(f,"%f\t",ch1IntensityStdDev);
	fprintf(f,"%f\t",ch1BackgroundStdDev);
	fprintf(f,"%f\t",ch1Diameter);
	fprintf(f,"%d\t",ch1Area);
	fprintf(f,"%f\t",ch1Footprint);
	fprintf(f,"%f\t",ch1Circularity);
	fprintf(f,"%f\t",ch1SpotUniformity);
	fprintf(f,"%f\t",ch1BkgUniformity);
	fprintf(f,"%f\t",ch1SignalNoiseRatio);
	fprintf(f,"%f\t",ch1Confidence);

	fprintf(f,"%f\t",ch2Intensity);
	fprintf(f,"%f\t",ch2Background);
	fprintf(f,"%f\t",ch2IntensityStdDev);
	fprintf(f,"%f\t",ch2BackgroundStdDev);
	fprintf(f,"%f\t",ch2Diameter);
	fprintf(f,"%d\t",ch2Area);
	fprintf(f,"%f\t",ch2Footprint);
	fprintf(f,"%f\t",ch2Circularity);
	fprintf(f,"%f\t",ch2SpotUniformity);
	fprintf(f,"%f\t",ch2BkgUniformity);
	fprintf(f,"%f\t",ch2SignalNoiseRatio);
	fprintf(f,"%f\t",ch2Confidence);

	fprintf(f,"%d\t",ch1DiameterFilter);
	fprintf(f,"%d\t",ch1AreaFilter);
	fprintf(f,"%d\t",ch1FootprintFilter);
	fprintf(f,"%d\t",ch1CircularityFilter);
	fprintf(f,"%d\t",ch1SpotUniformityFilter);
	fprintf(f,"%d\t",ch1BkgUniformityFilter);
	fprintf(f,"%d\t",ch1SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch1ReplicateFilter);

	fprintf(f,"%d\t",ch2DiameterFilter);
	fprintf(f,"%d\t",ch2AreaFilter);
	fprintf(f,"%d\t",ch2FootprintFilter);
	fprintf(f,"%d\t",ch2CircularityFilter);
	fprintf(f,"%d\t",ch2SpotUniformityFilter);
	fprintf(f,"%d\t",ch2BkgUniformityFilter);
	fprintf(f,"%d\t",ch2SignalNoiseRatioFilter);
	fprintf(f,"%d\t",ch2ReplicateFilter);

	fprintf(f,"%d\n",IgnoreFilter);

	return(IgnoreFilter);
};

void *ScanArraySpot::ValueOf( int code){
	switch(code){
		case CH1INTENSITY: return(&ch1Intensity);
   	case CH2INTENSITY: return(&ch2Intensity);
   	case CH1BACKGROUND: return(&ch1Background);
   	case CH2BACKGROUND: return(&ch2Background);
   	case CH1INTENSITYSTDEV: return(&ch1IntensityStdDev);
   	case CH2INTENSITYSTDEV: return(&ch2IntensityStdDev);
	   case CH1BACKGROUNDSTDEV: return(&ch1BackgroundStdDev);
   	case CH2BACKGROUNDSTDEV: return(&ch2BackgroundStdDev);
	   case CH1DIAMETER: return(&ch1Diameter);
	   case CH2DIAMETER: return(&ch2Diameter);
	   case CH1AREA: return(&ch1Area);
	   case CH2AREA: return(&ch2Area);
	   case CH1FOOTPRINT: return(&ch1Footprint);
	   case CH2FOOTPRINT: return(&ch2Footprint);
	   case CH1CIRCULARITY: return(&ch1Circularity);
	   case CH2CIRCULARITY: return(&ch2Circularity);
	   case CH1SPOTUNIFORMITY: return(&ch1SpotUniformity);
	   case CH2SPOTUNIFORMITY: return(&ch2SpotUniformity);
	   case CH1BKGUNIFORMITY: return(&ch1BkgUniformity);
	   case CH2BKGUNIFORMITY: return(&ch2BkgUniformity);
	   case CH1SIGNALNOISERATIO: return(&ch1SignalNoiseRatio);
	   case CH2SIGNALNOISERATIO: return(&ch2SignalNoiseRatio);
	   case CH1CONFIDENCE: return(&ch1Confidence);
	   case CH2CONFIDENCE: return(&ch2Confidence);
   }
	return(NULL);
}


ScanArraySpot::ScanArraySpot(FILE *f, int num){
	Read(f, num);
};

int ScanArraySpot::SummarizeFilters(){
   IgnoreFilter+=ch1DiameterFilter;
   IgnoreFilter+=ch1AreaFilter;
   IgnoreFilter+=ch1FootprintFilter;
   IgnoreFilter+=ch1CircularityFilter;
   IgnoreFilter+=ch1SpotUniformityFilter;
   IgnoreFilter+=ch1BkgUniformityFilter;
   IgnoreFilter+=ch1SignalNoiseRatioFilter;
   IgnoreFilter+=ch1ReplicateFilter;
   IgnoreFilter+=ch2DiameterFilter;
   IgnoreFilter+=ch2AreaFilter;
   IgnoreFilter+=ch2FootprintFilter;
   IgnoreFilter+=ch2CircularityFilter;
   IgnoreFilter+=ch2SpotUniformityFilter;
   IgnoreFilter+=ch2BkgUniformityFilter;
   IgnoreFilter+=ch2SignalNoiseRatioFilter;
   IgnoreFilter+=ch2ReplicateFilter;
	return(IgnoreFilter);

}

//////////////////////////////////////////////////////////////////////////////

ScanArraySlideData::ScanArraySlideData(){
	num_spots=0;
}

ScanArraySlideData::~ScanArraySlideData(){
int i;
	for(i=0;i<num_spots;i++){
   	if(spot[i]) delete spot[i];
   }
   if(header) delete header;
   if(protocol) delete protocol;
   if(tolerance) delete tolerance;
   if(image) delete image;
}

int ScanArraySlideData::Find(FILE *f){
char s[1000], *ss;
int flag;
	s[0]=0;
	flag=num_data_columns=0;
	while(flag<2){
   	if(!fgets(s,1000,f)){
      	rewind(f);
         flag++;
      }
      if(strstr(s,"Begin Data")) {
         fgets(s,1000,f);//Read and analyze the legend line
			ss=s;
  	   	while(1){
	      	ss=strchr(ss+1,'\t');
	         num_data_columns++;
				if(!ss) break;
	      }
      	return(flag);
      }
	}
   return(flag);
};

int ScanArraySlideData::Read(FILE *f){
int i;
ScanArraySpot *spt;
	if(!num_spots){
   	protocol = new ProtocolInfo();
      protocol->Find(f);
      protocol->Read(f);
   	num_spots=protocol->TotalSpots;
   }
	Find(f);
   //num_spots=10000;
   spot = new ScanArraySpot *[num_spots];// calloc(num_spots,sizeof(ScanArraySpot *));
   for(i=0;i<num_spots;i++){
   	spt= new ScanArraySpot(f, num_data_columns);
      spot[i]=spt;
   }
	return(1);
};

int ScanArraySlideData::Write(FILE *f){
int i;
	header->Write(f);
   protocol->Write(f);
   tolerance->Write(f);
   image->Write(f);
	fprintf(f,"Begin Data\n");
   fprintf(f,"Number\tArray Row\tArray Column\tRow\tColumn\tName\tX Location\tY Location\tch1 Intensity\tch1 Background\tch1 Intensity Std Dev\tch1 Background Std Dev\tch1 Diameter\tch1 Area\tch1 Footprint\tch1 Circularity\tch1 Spot Uniformity\tch1 Bkg. Uniformity\tch1 Signal Noise Ratio\tch1 Confidence\tch2 Intensity\tch2 Background\tch2 Intensity Std Dev\tch2 Background Std Dev\tch2 Diameter\tch2 Area\tch2 Footprint\tch2 Circularity\tch2 Spot Uniformity\tch2 Bkg. Uniformity\tch2 Signal Noise Ratio\tch2 Confidence\tch1 Diameter Filter\tch1 Area Filter\tch1 Footprint Filter\tch1 Circularity Filter\tch1 Spot Uniformity Filter\tch1 Bkg. Uniformity Filter\tch1 Signal Noise Ratio Filter\tch1 Replicate Filter\tch2 Diameter Filter\tch2 Area Filter\tch2 Footprint Filter\tch2 Circularity Filter\tch2 Spot Uniformity Filter\tch2 Bkg. Uniformity Filter\tch2 Signal Noise Ratio Filter\tch2 Replicate Filter\tIgnore Filter\n");

   for(i=0;i<num_spots;i++){
      spot[i]->Write(f);
   }
	fprintf(f,"End Data\n");
	return(1);
}

int ScanArraySlideData::WriteShort(FILE *f){
int i;
   fprintf(f,"Number\tArrayRow\tArrayColumn\tRow\tColumn\tName\tXLocation\tYLocation\tch1Intensity\tch1Background\tch1IntensityStdDev\tch1BackgroundStdDev\tch1Diameter\tch1Area\tch1Footprint\tch1Circularity\tch1SpotUniformity\tch1Bkg.Uniformity\tch1SignalNoiseRatio\tch1Confidence\tch2Intensity\tch2Background\tch2IntensityStdDev\tch2BackgroundStdDev\tch2Diameter\tch2Area\tch2Footprint\tch2Circularity\tch2SpotUniformity\tch2Bkg.Uniformity\tch2SignalNoiseRatio\tch2Confidence\tch1DiameterFilter\tch1AreaFilter\tch1FootprintFilter\tch1CircularityFilter\tch1SpotUniformityFilter\tch1Bkg.UniformityFilter\tch1SignalNoiseRatioFilter\tch1ReplicateFilter\tch2DiameterFilter\tch2AreaFilter\tch2FootprintFilter\tch2CircularityFilter\tch2SpotUniformityFilter\tch2Bkg.UniformityFilter\tch2SignalNoiseRatioFilter\tch2ReplicateFilter\tIgnoreFilter\n");
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f);
   }
	return(1);
}

int ScanArraySlideData::WriteShort(FILE *f, char *extraname, float *extravalue){
int i;
   fprintf(f,"Number\tArrayRow\tArrayColumn\tRow\tColumn\tName\tXLocation\tYLocation\t%s\tch1Intensity\tch1Background\tch1IntensityStdDev\tch1BackgroundStdDev\tch1Diameter\tch1Area\tch1Footprint\tch1Circularity\tch1SpotUniformity\tch1Bkg.Uniformity\tch1SignalNoiseRatio\tch1Confidence\tch2Intensity\tch2Background\tch2IntensityStdDev\tch2BackgroundStdDev\tch2Diameter\tch2Area\tch2Footprint\tch2Circularity\tch2SpotUniformity\tch2Bkg.Uniformity\tch2SignalNoiseRatio\tch2Confidence\tch1DiameterFilter\tch1AreaFilter\tch1FootprintFilter\tch1CircularityFilter\tch1SpotUniformityFilter\tch1Bkg.UniformityFilter\tch1SignalNoiseRatioFilter\tch1ReplicateFilter\tch2DiameterFilter\tch2AreaFilter\tch2FootprintFilter\tch2CircularityFilter\tch2SpotUniformityFilter\tch2Bkg.UniformityFilter\tch2SignalNoiseRatioFilter\tch2ReplicateFilter\tIgnoreFilter\n", extraname);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f,extravalue[i]);
   }
	return(1);
}

int ScanArraySlideData::WriteShort(FILE *f, char *extraname1, char *extraname2, char *extraname3, float *extravalue1, float *extravalue2, float *extravalue3){
int i;
   fprintf(f,"Number\tArrayRow\tArrayColumn\tRow\tColumn\tName\tXLocation\tYLocation\t%s\t%s\t%s\tch1Intensity\tch1Background\tch1IntensityStdDev\tch1BackgroundStdDev\tch1Diameter\tch1Area\tch1Footprint\tch1Circularity\tch1SpotUniformity\tch1Bkg.Uniformity\tch1SignalNoiseRatio\tch1Confidence\tch2Intensity\tch2Background\tch2IntensityStdDev\tch2BackgroundStdDev\tch2Diameter\tch2Area\tch2Footprint\tch2Circularity\tch2SpotUniformity\tch2Bkg.Uniformity\tch2SignalNoiseRatio\tch2Confidence\tch1DiameterFilter\tch1AreaFilter\tch1FootprintFilter\tch1CircularityFilter\tch1SpotUniformityFilter\tch1Bkg.UniformityFilter\tch1SignalNoiseRatioFilter\tch1ReplicateFilter\tch2DiameterFilter\tch2AreaFilter\tch2FootprintFilter\tch2CircularityFilter\tch2SpotUniformityFilter\tch2Bkg.UniformityFilter\tch2SignalNoiseRatioFilter\tch2ReplicateFilter\tIgnoreFilter\n", extraname1, extraname2, extraname3);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f,extravalue1[i],extravalue2[i], extravalue3[i]);
   }
	return(1);
}

int ScanArraySlideData::ReturnValuesFrom( ArrayStat *stat, int code){
int i,j,k;
int offset;
	offset=stat->offset;
	if((!stat->subarray)&&(offset)) return(0);
   switch (code){
	   case CH1INTENSITY:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->ch1Intensity=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case CH2INTENSITY:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->ch2Intensity=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case CH1BACKGROUND:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->ch1Background=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case CH2BACKGROUND:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->ch2Background=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   };
	return(code);
}

/////////////////////////////////////////////////////////////////////////////

AtlasArraySpot::AtlasArraySpot(){
	IgnoreFilter=1;
};

AtlasArraySpot::~AtlasArraySpot(){
	if(Genecode) delete(Genecode);
   if(ProteinGene) delete(ProteinGene);
};


int AtlasArraySpot::Read(FILE *f){
int i;
char s[1000],s1[100],s2[100];
char *ss;
//	i=0;
   fgets(s,1000,f);
   sscanf(s,"%s%d%d%d%d%d%d",s1,&Intensity_1,&Background_1,&AdjIntensity_1,&Intensity_2,&Background_2,&AdjIntensity_2);
	ss=s;
	for(i=0;i<9;i++) ss=(strchr(ss,'\t'))+1;
	Genecode=(char *)calloc(strlen(s1)+1, sizeof(char));
   strcpy(Genecode,s1);
   ProteinGene=(char *)calloc(strlen(ss)+1, sizeof(char));
   strcpy(ProteinGene,ss);
   if(AdjIntensity_1) Ratio=((float)AdjIntensity_2)/AdjIntensity_1;
   Difference=AdjIntensity_2-AdjIntensity_1;
	/*Genecode=FGetTabS(f); i++;
	i+=fscanf(f,"%d\t",&Intensity_1);
	i+=fscanf(f,"%d\t",&Background_1);
	i+=fscanf(f,"%d\t",&AdjIntensity_1);
	i+=fscanf(f,"%d\t",&Intensity_2);
	i+=fscanf(f,"%d\t",&Background_2);
	i+=fscanf(f,"%d\t",&AdjIntensity_2);
	i+=fscanf(f,"%f\t",&Ratio);
	i+=fscanf(f,"%d\t",&Difference);
	ProteinGene=FGetTabS(f); i++;*/
   IgnoreFilter=0;
   return(IgnoreFilter);
};

AtlasArraySpot::Write(FILE *f){
	fprintf(f,"%s\t",Genecode);
	fprintf(f,"%d\t",Intensity_1);
	fprintf(f,"%d\t",Background_1);
	fprintf(f,"%d\t",AdjIntensity_1);
	fprintf(f,"%d\t",Intensity_2);
	fprintf(f,"%d\t",Background_2);
	fprintf(f,"%d\t",AdjIntensity_2);
	fprintf(f,"%f\t",Ratio);
	fprintf(f,"%s\t",Difference);
	fprintf(f,"%s\n",ProteinGene);
	return(IgnoreFilter);
};

AtlasArraySpot::Write(FILE *f, float extravalue){
	fprintf(f,"%s\t",Genecode);
	fprintf(f,"%d\t",Intensity_1);
	fprintf(f,"%d\t",Background_1);
	fprintf(f,"%d\t",AdjIntensity_1);
	fprintf(f,"%d\t",Intensity_2);
	fprintf(f,"%d\t",Background_2);
	fprintf(f,"%d\t",AdjIntensity_2);
	fprintf(f,"%f\t",Ratio);
	fprintf(f,"%s\t",Difference);

   fprintf(f,"%f\t", extravalue);

	fprintf(f,"%s\n",ProteinGene);
	return(IgnoreFilter);
}


void *AtlasArraySpot::ValueOf( int code){
	switch(code){
   	case GENCODE: return(&Genecode);
		case INTENSITY_1: return(&Intensity_1);
   	case BACKGROUND_1: return(&Background_1);
		case ADJINTENSITY_1: return(&AdjIntensity_1);
		case INTENSITY_2: return(&Intensity_2);
   	case BACKGROUND_2: return(&Background_2);
		case ADJINTENSITY_2: return(&AdjIntensity_2);
   	case RATIO: return(&Ratio);
   	case DIFFERENCE: return(&Difference);
   	case PROTEINGENE: return(&ProteinGene);
   }
	return(NULL);
}


AtlasArraySpot::AtlasArraySpot(FILE *f){
	Read(f);
};


/////////////////////////////////////////////////////////////////////////////

AtlasArraySlideData::AtlasArraySlideData(){
	num_spots=0;
}

AtlasArraySlideData::~AtlasArraySlideData(){
int i;
	for(i=0;i<num_spots;i++){
   	if(spot[i]) delete spot[i];
   }
   if(header) delete header;
}

int AtlasArraySlideData::Find(FILE *f){
	rewind(f);
   return(1);
};

int AtlasArraySlideData::Read(FILE *f){
int i;
char s[1000];
AtlasArraySpot *spt;
	Find(f);
	if(!num_spots){
	   fgets(s,1000,f);
      while(fgets(s,1000,f)){
      	if(isalnum(s[0])) num_spots++;
         s[0]=0;
      }
		Find(f);
   }
	fgets(s,1000,f);
   header=(char *)calloc(strlen(s)+1,sizeof(char));
   strcpy(header,s);
	fgets(s,1000,f);

   spot = new AtlasArraySpot *[num_spots];// calloc(num_spots,sizeof(ScanArraySpot *));
   for(i=0;i<num_spots;i++){
   	spt= new AtlasArraySpot(f);
      spot[i]=spt;
   }
	return(1);
};

int AtlasArraySlideData::Read(FILE *f, int layout){
int i;
char s[1000];
AtlasArraySpot *spt;
	Find(f);
	if(layout==PLASTIC_RAT4K){
   	num_data_rows=8;
      num_data_columns=3;
      num_array_rows=16;
      num_array_columns=24;
   }
	if(!num_spots){
	   fgets(s,1000,f);
      while(fgets(s,1000,f)){
      	if(isalnum(s[0])) num_spots++;
         s[0]=0;
      }
		Find(f);
   }
	fgets(s,1000,f);
   header=(char *)calloc(strlen(s)+1,sizeof(char));
   strcpy(header,s);
	fgets(s,1000,f);

   spot = new AtlasArraySpot *[num_spots];// calloc(num_spots,sizeof(ScanArraySpot *));
   for(i=0;i<num_spots;i++){
   	spt= new AtlasArraySpot(f);
      spot[i]=spt;
   }
	return(1);
};

int AtlasArraySlideData::Write(FILE *f){
int i;
	fprintf(f,"%s",header);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f);
   }
	return(1);
}


int AtlasArraySlideData::Write(FILE *f, char *extraname, float *extravalue){
int i;
   fprintf(f,"Gene code\tIntensity_1\tBackground_1\tAdj.Intensity_1\tIntensity_2\tBackground_2\tAdj.Intensity_2\tRatio\tDifference\t%s\tProtein/gene\n",extraname);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f,extravalue[i]);
   }
	return(1);
}

int AtlasArraySlideData::ReturnValuesFrom( ArrayStat *stat, int code){
int i,j,k;
int offset;
	offset=stat->offset;
	if((!stat->subarray)&&(offset)) return(0);
   switch (code){
	   case INTENSITY_1:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Intensity_1=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case INTENSITY_2:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Intensity_2=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
	   case ADJINTENSITY_1:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->AdjIntensity_1=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case ADJINTENSITY_2:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->AdjIntensity_2=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case BACKGROUND_1:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Background_1=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case BACKGROUND_2:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Background_2=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case RATIO:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Ratio=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case DIFFERENCE:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Difference=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   };
	return(code);
}

ArrayStat *AtlasArraySlideData::MakeStat(int layout, int code){
int i,row,column;
int p;
ArrayStat *stat;
	stat= new ArrayStat();
   switch(layout){
   	case PLASTIC_RAT4K:{
      	stat->code=PLASTIC_RAT4K;
         stat->subarray=384;
         stat->offset=0;
         stat->numRows=128;
         stat->numColumns=72;
			stat->matrix=(float **)calloc(stat->numRows, sizeof(float *));
			for(i=0;i<stat->numRows;i++) stat->matrix[i]=(float *)calloc(stat->numColumns,sizeof(float));
         for(i=0;i<num_spots;i++){
				if(spot[i]->IgnoreFilter) continue;
            row=(int)(spot[i]->Genecode[0]-'A')*8+(int)(spot[i]->Genecode[5]-'1');
            column=(int)(spot[i]->Genecode[1]-'0')*10 + (int)(spot[i]->Genecode[2]-'1');
            column=column*3;
				if(spot[i]->Genecode[3]=='c') column++;
				if(spot[i]->Genecode[3]=='e') column+=2;
				if(code==RATIO){
   	         stat->matrix[row][column]=*(float *)(spot[i]->ValueOf(code));
				}else{
	            p=*(int *)(spot[i]->ValueOf(code));
   	         stat->matrix[row][column]=(float)p;
            }
         }
      };break;
   }
   return(stat);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

GenericArraySpot::GenericArraySpot(){
	IgnoreFilter=1;
};

GenericArraySpot::~GenericArraySpot(){
	if(Name) delete(Name);
   if(ProteinGene) delete(ProteinGene);
};


int GenericArraySpot::Read(FILE *f){
int i;
char s[1000],s1[100];
char *ss;
   fgets(s,1000,f);
   sscanf(s,"%s%d%d%d%d",s1,&Intensity_1,&Background_1,&Intensity_2,&Background_2);
	ss=s;
	for(i=0;i<9;i++) ss=(strchr(ss,'\t'))+1;
	Name=(char *)calloc(strlen(s1)+1, sizeof(char));
   strcpy(Name,s1);
   if(Intensity_1) Ratio=((float)Intensity_2)/Intensity_1;
   Difference=Intensity_2-Intensity_1;
   IgnoreFilter=0;
   return(IgnoreFilter);
};

GenericArraySpot::Write(FILE *f){
	fprintf(f,"%s\t",Name);
	fprintf(f,"%d\t",Intensity_1);
	fprintf(f,"%d\t",Background_1);
	fprintf(f,"%d\t",Intensity_2);
	fprintf(f,"%d\t",Background_2);
	fprintf(f,"%f\t",Ratio);
	fprintf(f,"%s\t",Difference);
	return(IgnoreFilter);
};

GenericArraySpot::Write(FILE *f, float extravalue){
	fprintf(f,"%s\t",Name);
	fprintf(f,"%d\t",Intensity_1);
	fprintf(f,"%d\t",Background_1);
	fprintf(f,"%d\t",Intensity_2);
	fprintf(f,"%d\t",Background_2);
	fprintf(f,"%f\t",Ratio);
	fprintf(f,"%s\t",Difference);

   fprintf(f,"%f\t", extravalue);

	return(IgnoreFilter);
}


void *GenericArraySpot::ValueOf( int code){
	switch(code){
   	case GENCODE: return(&Name);
		case INTENSITY_1: return(&Intensity_1);
   	case BACKGROUND_1: return(&Background_1);
		case INTENSITY_2: return(&Intensity_2);
   	case BACKGROUND_2: return(&Background_2);
   	case RATIO: return(&Ratio);
   	case DIFFERENCE: return(&Difference);
   }
	return(NULL);
}


GenericArraySpot::GenericArraySpot(FILE *f){
	Read(f);
};


/////////////////////////////////////////////////////////////////////////////

GenericArraySlideData::GenericArraySlideData(){
	num_spots=0;
}

GenericArraySlideData::~GenericArraySlideData(){
int i;
	for(i=0;i<num_spots;i++){
   	if(spot[i]) delete spot[i];
   }
   if(header) delete header;
}

int GenericArraySlideData::Find(FILE *f){
	rewind(f);
   return(1);
};

int GenericArraySlideData::ReadLayout(FILE *f){
int i;
char s[1000];
//char *ss, *sss;
//GenericArraySpot *spt;
//the paprameter file can be used to read the pingroup layout
	Find(f);
   i=0;
	while(fgets(s,1000,f)){
   	*strchr(s,'\n')=0;
	   if(strstr(s,"arrayrows=")){ //rows of pigroups
			//default 1
	      num_array_rows=1;
	   	sscanf(strchr(s,'=')+1,"%d",&(num_array_rows));
         i++;
         continue;
	   }
	   if(strstr(s,"arraycolumns=")){ //columns of pingroups
			//default 0
	      num_array_columns=1;
	   	sscanf(strchr(s,'=')+1,"%d",&(num_array_columns));
         i++;
         continue;
	   }
	   if(strstr(s,"rows=")){ //rows of spots per pigroup
			//default 0
	      num_data_rows=0;
	   	sscanf(strchr(s,'=')+1,"%d",&(num_data_rows));
         i++;
         continue;
	   }
	   if(strstr(s,"columns=")){ //columns of spots per pingroup
			//default 0
	      num_data_columns=0;
	   	sscanf(strchr(s,'=')+1,"%d",&(num_data_columns));
         i++;
         continue;
	   }
	   if(strstr(s,"totalspots=")){ //total number of spots
			//default 0
	      num_spots=0;
	   	sscanf(strchr(s,'=')+1,"%d",&(num_spots));
         i++;
         continue;
	   }
	   if(strstr(s,"ch1intensity=")){ //which column is ch1 intensity?
			//default 1
	      ch1intensity=1;
	   	sscanf(strchr(s,'=')+1,"%d",&(ch1intensity));
         continue;
	   }
	   if(strstr(s,"ch1background=")){ //which column is ch1 background?
			//default 2
	      ch1background=2;
	   	sscanf(strchr(s,'=')+1,"%d",&(ch1background));
         continue;
	   }
	   if(strstr(s,"ch2intensity=")){ //which column is ch2 intensity?
			//default 3
	      ch2intensity=3;
	   	sscanf(strchr(s,'=')+1,"%d",&(ch1intensity));
         continue;
	   }
	   if(strstr(s,"ch2background=")){ //which column is ch2 background?
			//default 4
	      ch2background=4;
	   	sscanf(strchr(s,'=')+1,"%d",&(ch2background));
         continue;
	   }
   }
	return(i);
}

int GenericArraySlideData::Read(FILE *f){
int i;
char s[1000];
//char *ss, *sss;
GenericArraySpot *spt;
	if(!num_spots){
	   fgets(s,1000,f);
      while(fgets(s,1000,f)){
      	if(isalnum(s[0])) num_spots++;
         s[0]=0;
      }
		Find(f);
   }
	//first line is presumed to be a header line with column names
	fgets(s,1000,f);
   header=(char *)calloc(strlen(s)+1,sizeof(char));
   strcpy(header,s);
	fgets(s,1000,f);

   spot = new GenericArraySpot *[num_spots];// calloc(num_spots,sizeof(ScanArraySpot *));
   for(i=0;i<num_spots;i++){
   	spt= new GenericArraySpot(f);
      spot[i]=spt;
   }
	return(1);
};

int GenericArraySlideData::Write(FILE *f){
int i;
	fprintf(f,"%s",header);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f);
   }
	return(1);
}

int GenericArraySlideData::Write(FILE *f, char *extraname, float *extravalue){
int i;
   fprintf(f,"Gene code\tIntensity_1\tBackground_1\t\tIntensity_2\tBackground_2\tRatio\tDifference\t%s\tProtein/gene\n",extraname);
   for(i=0;i<num_spots;i++){
      spot[i]->Write(f,extravalue[i]);
   }
	return(1);
}

int GenericArraySlideData::ReturnValuesFrom( ArrayStat *stat, int code){
int i,j,k;
int offset;
	offset=stat->offset;
	if((!stat->subarray)&&(offset)) return(0);
   switch (code){
	   case INTENSITY_1:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Intensity_1=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case INTENSITY_2:{
   		k=0;
  			for(i=0;i<stat->numRows;i++){
	   		for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Intensity_2=stat->matrix[i][j];
               k++;
	         }
   	   }
	   };break;
   	case BACKGROUND_1:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Background_1=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case BACKGROUND_2:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Background_2=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case RATIO:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Ratio=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   	case DIFFERENCE:{
   		k=0;
	  		for(i=0;i<stat->numRows;i++){
		   	for(j=0;j<stat->numColumns;j++){
					spot[offset+k]->Difference=stat->matrix[i][j];
               k++;
         	}
	      }
	   };break;
   };
	return(code);
}

ArrayStat *GenericArraySlideData::MakeStat(int field){
//codes: 0- intensity1 1-backgrnd1 2- intensity2 3- backgrnd2 4- ratio
int i,j,k;
ArrayStat *stat;
	stat= new ArrayStat();
	stat->matrix=(float **)calloc(num_data_rows, sizeof(float *));
	//for the whole slide, filled from a particular field
	stat->code=field;
	stat->name=NULL;
	stat->number=0;
	stat->Eigval=NULL;
	stat->Eigvec=NULL;
	stat->numRows=num_array_rows*num_data_rows;
	stat->numColumns=num_array_columns*num_data_columns;
	stat->matrix=(float **)calloc(stat->numRows, sizeof(float *));
	for(i=0;i<stat->numRows;i++) stat->matrix[i]=(float *)calloc(stat->numColumns,sizeof(float));
	k=0;
	for(i=0;i<stat->numRows;i++){
 		for(j=0;j<stat->numColumns;j++){
			stat->matrix[i][j]=*(float *)(spot[k++]->ValueOf(field));
      }
   }
   return(stat);
}

ArrayStat *GenericArraySlideData::MakeStat(int field, int arraynumber){
//codes: 0- intensity1 1-backgrnd1 2- intensity2 3- backgrnd2 4- ratio
int i,j,k;
ArrayStat *stat;
	stat= new ArrayStat();
	stat->matrix=(float **)calloc(num_data_rows, sizeof(float *));
	stat->code=field;
	stat->name=NULL;
	stat->number=0;
   stat->subarray=1;
	stat->Eigval=NULL;
	stat->Eigvec=NULL;
	stat->numColumns=num_data_columns*num_array_columns;
	stat->numRows=(num_spots/3) / num_data_columns;

	stat->matrix=(float **)calloc(stat->numRows, sizeof(float *));
	for(i=0;i<stat->numRows;i++) stat->matrix[i]=(float *)calloc(stat->numColumns,sizeof(float));

   stat->offset=(num_spots/3)*arraynumber;
	k=0;
	for(i=0;i<stat->numRows;i++){
  		for(j=0;j<stat->numColumns;j++){
			stat->matrix[i][j]=*(float *)(spot[stat->offset+k]->ValueOf(field));
         k++;
      }
   }
	return(stat);
}

ArrayStat *GenericArraySlideData::MakeStat(int field, int subarray_row, int subarray_col){
int i,j,k;
ArrayStat *stat;
	stat= new ArrayStat();
	stat->matrix=(float **)calloc(num_data_rows, sizeof(float *));
	stat->code=field;
	stat->name=NULL;
	stat->number=0;
   stat->subarray=1;
	stat->Eigval=NULL;
	stat->Eigvec=NULL;
	stat->numRows=num_data_rows;
	stat->numColumns=num_data_columns;
   
	stat->matrix=(float **)calloc(stat->numRows, sizeof(float *));
	for(i=0;i<stat->numRows;i++) stat->matrix[i]=(float *)calloc(stat->numColumns,sizeof(float));

   k=num_data_columns*num_data_rows;
   stat->offset=k*num_array_columns*subarray_row+k*subarray_col;
	k=0;
	for(i=0;i<stat->numRows;i++){
  		for(j=0;j<stat->numColumns;j++){
			stat->matrix[i][j]=*(float *)(spot[stat->offset+k]->ValueOf(field));
         k++;
      }
   }
	return(stat);
}

/////////////////////////////////////////////////////////////////////////////


ArrayStat::ArrayStat(){
	code=0;
	subarray=0;
	offset=0;
	Mean=0;
   Median=0;
   Stdev=0;
};

ArrayStat::ArrayStat(ArrayStat *stat){
int i;
	code=stat->code;
	subarray=stat->subarray;
	offset=stat->offset;
	Mean=stat->Mean;
   Median=stat->Median;
   Stdev=stat->Stdev;
   numRows=stat->numRows;
   numColumns=stat->numColumns;
   matrix=(float **)malloc(numRows*sizeof(float *));
   for(i=0;i<numColumns;i++){
   	matrix[i]=(float *)malloc(numColumns*sizeof(float));
      memcpy(matrix[i],stat->matrix[i],numColumns*sizeof(float));
   }
};


ArrayStat::ArrayStat(ScanArraySlideData *slide){
	//for the whole slide, empty
	MatrixObject(slide->protocol->ArrayRows*slide->protocol->Rows, slide->protocol->ArrayColumns*slide->protocol->Columns);
};

ArrayStat::ArrayStat(ScanArraySlideData *slide, int field){
	//for the whole slide, filled from a particular field
int i,j,k;
	code=field;
	name=NULL;
	offset=0;
	number=0;
   subarray=0;
   precision=0;
	Eigval=NULL;
	Eigvec=NULL;
	numRows=slide->protocol->ArrayRows*slide->protocol->Rows;
	numColumns=slide->protocol->ArrayColumns*slide->protocol->Columns;
	matrix=(float **)calloc(numRows, sizeof(float *));
	for(i=0;i<numRows;i++) matrix[i]=(float *)calloc(numColumns,sizeof(float));
	k=0;
	for(i=0;i<numRows;i++){
 		for(j=0;j<numColumns;j++){
			matrix[i][j]=*(float *)(slide->spot[k++]->ValueOf(code));
      }
   }
};

ArrayStat::ArrayStat(ScanArraySlideData *slide, int field, int arraynumber){
int i,j,k;
	code=field;
	name=NULL;
	number=0;
   subarray=1;
	offset=0;
   precision=0;
	Eigval=NULL;
	Eigvec=NULL;
	numColumns=slide->protocol->Columns*slide->protocol->ArrayColumns;
	numRows=(slide->protocol->TotalSpots/3) / numColumns;
	matrix=(float **)calloc(numRows, sizeof(float *));
	for(i=0;i<numRows;i++) matrix[i]=(float *)calloc(numColumns,sizeof(float));

   offset=(slide->protocol->TotalSpots/3)*arraynumber;
	k=0;
	for(i=0;i<numRows;i++){
  		for(j=0;j<numColumns;j++){
			matrix[i][j]=*(float *)(slide->spot[offset+k]->ValueOf(code));
         k++;
      }
   }
};


ArrayStat::ArrayStat(ScanArraySlideData *slide, int field, int subarray_row, int subarray_col){
int i,j,k;
	code=field;
	name=NULL;
	number=0;
   subarray=1;
	Eigval=NULL;
	Eigvec=NULL;
	numRows=slide->protocol->Rows;
	numColumns=slide->protocol->Columns;
	matrix=(float **)calloc(numRows, sizeof(float *));
	for(i=0;i<numRows;i++) matrix[i]=(float *)calloc(numColumns,sizeof(float));

   k=slide->protocol->Columns*slide->protocol->Rows;
   offset=k*slide->protocol->ArrayColumns*subarray_row+k*subarray_col;
	k=0;
	for(i=0;i<numRows;i++){
  		for(j=0;j<numColumns;j++){
			matrix[i][j]=*(float *)(slide->spot[offset+k]->ValueOf(code));
         k++;
      }
   }
};

ArrayStat::~ArrayStat(){
};


float ArrayStat::get_mean(){
int i,j;
	Mean=0;
	for(i=0;i<numRows;i++){
   	for(j=0;j<numColumns;j++){
      	Mean+=matrix[i][j];
      }
   }
   Mean=Mean/(numRows*numColumns);
   return(Mean);
}

float ArrayStat::get_median(){
int i,j,k;
float *a;
	k=0;
	a=(float *)calloc(numRows*numColumns,sizeof(float));
	for(i=0;i<numRows;i++){
   	for(j=0;j<numColumns;j++){
      	a[k++]=matrix[i][j];
      }
   }
   qsort((void *)a, k, sizeof(float),float_sorter);
   Median=a[k/2];
   free(a);
   return(Median);
}

float ArrayStat::get_MAD(){
int i,j,k;
float *a;
	if(!Median) get_median();
   if(!Median) return(0);
	for(i=0;i<numRows;i++){
   	for(j=0;j<numColumns;j++){
      	a[k++]=matrix[i][j]-Median;
      }
   }
   qsort((void *)a, k, sizeof(float),float_sorter);
   Median=a[k/2];
   free(a);
   return(MAD);
}

float ArrayStat::get_TMAD(int percentile){
int i,j,k;
float x;
float *a;
	if(!Mean) get_trimmed_mean(percentile);
   if(!Mean) return(0);
	a=(float *)calloc(numRows*numColumns,sizeof(float));
	for(i=0;i<numRows;i++){
   	for(j=0;j<numColumns;j++){
      	a[k++]=matrix[i][j]-Mean;
      }
   }
   qsort((void *)a, k, sizeof(float),float_sorter);
   Median=a[k/2];
	k=(numRows*numColumns*percentile)/100;
	x=0;
   for(i=k;i<numRows*numColumns-k;i++){
   	x+=a[i];
   }
   Mean=x/(numRows*numColumns-(k*2));
   free(a);
   return(TMAD);
}

float ArrayStat::get_trimmed_mean(int percentile){
int i,j,k;
float x;
float *a;
	k=0;
	a=(float *)calloc(numRows*numColumns,sizeof(float));
	for(i=0;i<numRows;i++){
   	for(j=0;j<numColumns;j++){
      	a[k++]=matrix[i][j];
      }
   }
   qsort((void *)a, k, sizeof(float),float_sorter);
   Median=a[k/2];
	k=(numRows*numColumns*percentile)/100;
	x=0;
   for(i=k;i<numRows*numColumns-k;i++){
   	x+=a[i];
   }
   Mean=x/(numRows*numColumns-(k*2));
   free(a);
   return(Mean);
}

float ArrayStat::get_stdev(){
int i,j;
	Stdev=0;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			Stdev+=(Mean - matrix[i][j])*(Mean - matrix[i][j]);
      }
   }
	Stdev=sqrt(Stdev/((numRows*numColumns)-1));
   return(Stdev);
}

void ArrayStat::ScaleToMean(){
int i,j;
	if(!Mean) get_mean();
   if(!Mean) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/Mean;
      }
   }
}

void ArrayStat::ScaleToStdev(){
int i,j;
	if(!Stdev) get_stdev();
   if(!Stdev) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/Stdev;
      }
   }
}
void ArrayStat::ZScore(){
int i,j;
	if(!Mean) get_mean();
   if(!Stdev) get_stdev();
   if(!Stdev) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=(matrix[i][j]-Mean)/Stdev;
      }
   }
}

void ArrayStat::ScaleToMedian(){
int i,j;
	if(!Median) get_median();
   if(!Median) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/Median;
      }
   }
}

void ArrayStat::ScaleToMAD(){
int i,j;
	if(!MAD) get_MAD();
   if(!MAD) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/MAD;
      }
   }
}

void ArrayStat::ScaleToTMAD(int percentile){
int i,j;
	if(!TMAD) get_TMAD(percentile);
   if(!TMAD) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/TMAD;
      }
   }
}

void ArrayStat::GetLog2(){
int i,j;
	if(!Median) get_mean();
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=logl(matrix[i][j])/logl(2);
      }
   }
}

void ArrayStat::ScaleTo(float factor){
int i,j;
	if(factor==0) return;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			matrix[i][j]=matrix[i][j]/factor;
      }
   }
}

void ArrayStat::shuffle_by_xy(int times){
int i,x,y,x1,y1;
float b;
	for(i=0;i<times;i++){
   	x=rand()%numRows;
      y=rand()%numColumns;
      x1=rand()%numRows;
      y1=rand()%numColumns;
      b=matrix[x][y];
      matrix[x][y]=matrix[x1][y1];
      matrix[x1][y1]=b;
   }
}

int ArrayStat::mean_equal_to(ArrayStat *other, int bon, float *d){
int j;
float t,s2;
	Mean=get_trimmed_mean(5);
   //Mean=get_mean();
	other->Mean=other->get_trimmed_mean(5);
   //other->Mean=other->get_mean();
   Stdev=get_stdev();
   other->Stdev=other->get_stdev();
	s2=Stdev*(numRows * numColumns - 1)+other->Stdev*(other->numRows * other->numColumns - 1);
	s2=s2/((numRows * numColumns)+(other->numRows * other->numColumns)-2);
   t=fabs((Mean-other->Mean)/sqrt(s2/(numRows * numColumns)+s2/(other->numRows * other->numColumns)));
	if(!bon) bon=1;
   j=(numRows * numColumns - 1)+(other->numRows * other->numColumns - 1);
   *d=(F_table(1,j)/bon)-t;
   if( t >(F_table(1,j)/bon)) {
   	return(0);
   }else{
   	return(1);
   }
}

int ArrayStat::grid_is_uniform( int grid ){
int i,j,k,l,xstep,ystep;
float a;
ArrayStat *part;

	xstep=numRows/grid; ystep=numColumns/grid;
   part=new ArrayStat();
   part->matrix=(float **)calloc(xstep,sizeof(float *));
	part->numRows=xstep; part->numColumns=ystep;
   for(i=0;i<xstep;i++) part->matrix[i]=(float *)calloc(ystep,sizeof(float));

   for(i=0;i<xstep*grid;i+=xstep){
   	for(j=0;j<ystep*grid;j+=ystep){
      	for(k=0;k<xstep;k++){
         	for(l=0;l<ystep;l++){
            	part->matrix[k][l]=matrix[i+k][j+l];
            }
         }
			/*with positive 'a' means are equal*/
         if(!mean_equal_to(part, 1, &a)) {
				delete(part);
         	return(0);
         }
      }
   }
   return(1);
}

int ArrayStat::grid_is_bonferroni_uniform( int grid ){
int i,j,k,l,xstep,ystep;
float a;
ArrayStat *part;

	xstep=numRows/grid; ystep=numColumns/grid;
   part=new ArrayStat();
   part->matrix=(float **)calloc(xstep,sizeof(float *));
	part->numRows=xstep; part->numColumns=ystep;
   for(i=0;i<xstep;i++) part->matrix[i]=(float *)calloc(ystep,sizeof(float));

   for(i=0;i<xstep*grid;i+=xstep){
   	for(j=0;j<ystep*grid;j+=ystep){
      	for(k=0;k<xstep;k++){
         	for(l=0;l<ystep;l++){
            	part->matrix[k][l]=matrix[i+k][j+l];
            }
         }
			/*with positive 'a' means are equal*/
         if(!mean_equal_to(part, grid*grid, &a)) {
				delete(part);
         	return(0);
         }
      }
   }
   return(1);
}


int ArrayStat::Write(FILE *f){
int i,j;
   for(i=0;i<numRows; i++){
   	for(j=0;j<numColumns;j++){
			fprintf(f,"%f\d",matrix[i][j]);
      }
      fprintf(f,"\n");
   }
   return(1);
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
VectorObject *GetSlideControls(int code, char *keyword, ScanArraySlideData *slide){
int i,j;
ScanArraySpot *spot;
VectorObject *controls;
	j=0;
	for(i=0;i<slide->protocol->TotalSpots;i++){
   	if(strstr(slide->spot[i]->Name, keyword)) j++;
   }
	controls=new VectorObject(j, FLOAT);
	j=0;
	for(i=0;i<slide->protocol->TotalSpots;i++){
		if(strstr(slide->spot[i]->Name, keyword)) {
      	spot=slide->spot[i];
      	controls->vec[j++]=*(float *)(spot->ValueOf(code));
      }
   }
   return(controls);
};

VectorObject *GetSlideControls(int code, char *keyword, GenericArraySlideData *slide, int TotalSpots){
int i,j;
GenericArraySpot *spot;
VectorObject *controls;
	j=0;
	for(i=0;i<TotalSpots;i++){
   	if(strstr(slide->spot[i]->Name, keyword)) j++;
   }
	controls=new VectorObject(j, FLOAT);
	j=0;
	for(i=0;i<TotalSpots;i++){
		if(strstr(slide->spot[i]->Name, keyword)) {
      	spot=slide->spot[i];
      	controls->vec[j++]=*(float *)(spot->ValueOf(code));
      }
   }
   return(controls);
};

VectorObject *GetSlideSubarrayControls(int code, char *keyword, ScanArraySlideData *slide, int ArrayRow, int ArrayColumn){
int i,j;
ScanArraySpot *spot;
VectorObject *controls;
	j=0;
	for(i=0;i<slide->protocol->TotalSpots;i++){
   	if(strstr(slide->spot[i]->Name, keyword)) j++;
   }
	controls=new VectorObject(j, FLOAT);
	j=0;
	for(i=0;i<slide->protocol->TotalSpots;i++){
     	spot=slide->spot[i];
		if(spot->IgnoreFilter) continue;
		if(spot->ArrayRow!=(ArrayRow+1)) continue;
      if(spot->ArrayColumn!=(ArrayColumn+1)) continue;
		if(strstr(slide->spot[i]->Name, keyword)) {
      	controls->vec[j++]=*(float *)(spot->ValueOf(code));
      }
   }
   controls->number=j;
   return(controls);
};

VectorObject *GetSlideSubarrayControls(int code, char *keyword, GenericArraySlideData *slide, int TotalSpots, int ArrayRow, int ArrayColumn){
int i,j;
GenericArraySpot *spot;
VectorObject *controls;
	j=0;
	for(i=0;i<TotalSpots;i++){
   	if(strstr(slide->spot[i]->Name, keyword)) j++;
   }
	controls=new VectorObject(j, FLOAT);
	j=0;
	for(i=0;i<TotalSpots;i++){
     	spot=slide->spot[i];
		if(spot->IgnoreFilter) continue;
		if(spot->ArrayRow!=(ArrayRow+1)) continue;
      if(spot->ArrayColumn!=(ArrayColumn+1)) continue;
		if(strstr(slide->spot[i]->Name, keyword)) {
      	controls->vec[j++]=*(float *)(spot->ValueOf(code));
      }
   }
   controls->number=j;
   return(controls);
};

int ld_sorter( const void *x, const void *y){
LowessData *xx,*yy;
	xx=*(LowessData **)x; yy=*(LowessData **)y;
	return((int)(xx->intensity - yy->intensity));
}


int ld_restorer(const void *x, const void *y){
LowessData *xx,*yy;
	xx=*(LowessData **)x; yy=*(LowessData **)y;
	return((int)(xx->number-yy->number));
}


int lowess_write(LowessData **ld, int num_spots, FILE *f){
int i;
	fprintf(f,"Log2(RG)\tLog2(R/G)\tLog2(R/G)_adjusted\tPredicted\n");
	for(i=0;i<num_spots;i++){
   	fprintf(f,"%f\t%f\t%f\t%f\n", ld[i]->intensity, ld[i]->lr, ld[i]->lr - ld[i]->logratio, ld[i]->logratio);
   }
   fprintf(f,"\n");
   return(1);
}

ArrayStat *Slidemean(ArrayStat *arr1, ArrayStat *arr2){
int i,j,k, row, column, flag, num_spots;
LowessData **ld;
float min, min_diff;//, middle;
float *buf, *buf1, *thisbuffer, *thatbuffer;

	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	buf1=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=log10(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf1[i]=buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
   	min=fabs(buf[i]-buf[i+1]);
      if(min==0) continue;
      if(min_diff>min) min_diff=min;
   }
	thisbuffer=buf;
   thatbuffer=buf1;
	while(1){
		thisbuffer[0]=(thisbuffer[0]+thisbuffer[1])/2;
	   for(i=1;i<num_spots-1;i++){
      	thisbuffer[i]=(thisbuffer[i-1]+thisbuffer[i]+thisbuffer[i+1])/3;
      }
      thisbuffer[i]=(thisbuffer[i]+thisbuffer[i-1])/2;
		flag=1;
	   for(i=0;i<num_spots;i++){
   		if(fabs(thisbuffer[i]-thatbuffer[i])>min_diff){
         	flag=0;
            break;
			}
      }
      if(flag){
      	break;
      }else{
		   memcpy(thatbuffer,thisbuffer,num_spots*sizeof(float));
      }
   }
//	middle=logl((arr1->get_trimmed_mean(10)+MINPOS)/(arr2->get_trimmed_mean(10)+MINPOS))/logl(2);
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	//qsort(ld, num_spots, sizeof(LowessData **), ld_restorer);
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio/* + middle*/;
   }
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
   free(buf1);
	return(arr1);
}

//smoothing with 3rd order local polynomial fit
ArrayStat *Polysmooth(ArrayStat *arr1, ArrayStat *arr2){
int i,j,k, row, column, num_spots;
LowessData **ld;
float min, min_diff;
float *buf, *thisbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thisbuffer[0]=(39*thisbuffer[0]+8*thisbuffer[1]-4*(thisbuffer[2]+thisbuffer[3]+thisbuffer[4])+thisbuffer[5]-2*thisbuffer[6])/42;
	thisbuffer[1]=(8*thisbuffer[0]+19*thisbuffer[1]+16*thisbuffer[2]+6*thisbuffer[3]-4*thisbuffer[4]-7*thisbuffer[5]+4*thisbuffer[7])/42;
	thisbuffer[2]=(-4*thisbuffer[0]+16*thisbuffer[1]+19*thisbuffer[2]+12*thisbuffer[3]+2*thisbuffer[4]-4*thisbuffer[5]+thisbuffer[6])/42;
   for(i=3;i<num_spots-3;i++){
     	thisbuffer[i]=(7*thisbuffer[i]+6*(thisbuffer[i+1]+thisbuffer[i-1])+3*(thisbuffer[i+2]+thisbuffer[i-2])-2*(thisbuffer[i+3]+thisbuffer[i-3]))/21;
   }
	thisbuffer[num_spots-3]=(thisbuffer[num_spots-7]-4*thisbuffer[num_spots-6]+2*thisbuffer[num_spots-5]+12*thisbuffer[num_spots-4]+19*thisbuffer[num_spots-3]+16*thisbuffer[num_spots-2]-4*thisbuffer[num_spots-1])/42;
	thisbuffer[num_spots-2]=(4*thisbuffer[num_spots-7]-7*thisbuffer[num_spots-6]-4*thisbuffer[num_spots-5]+6*thisbuffer[num_spots-4]+16*thisbuffer[num_spots-3]+19*thisbuffer[num_spots-2]+8*thisbuffer[num_spots-1])/42;
	thisbuffer[num_spots-1]=(-2*thisbuffer[num_spots-7]+4*thisbuffer[num_spots-6]+thisbuffer[num_spots-5]-4*thisbuffer[num_spots-4]-4*thisbuffer[num_spots-3]+8*thisbuffer[num_spots-2]+39*thisbuffer[num_spots-1])/42;
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
	return(arr1);
}

//smoothing with local linear fit
ArrayStat *Lsmooth(ArrayStat *arr1, ArrayStat *arr2){
int i,j,k, row, column, num_spots;
LowessData **ld;
float min, min_diff;
float *buf, *thisbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thisbuffer[0]=(3*thisbuffer[0]+2*thisbuffer[1]+thisbuffer[2]-thisbuffer[4])/5;
	thisbuffer[1]=(4*thisbuffer[0]+3*thisbuffer[1]+2*thisbuffer[2]+thisbuffer[3])/10;
   for(i=3;i<num_spots-2;i++){
     	thisbuffer[i]=(thisbuffer[i-2]+thisbuffer[i-1]+thisbuffer[i]+thisbuffer[i+1]+thisbuffer[i+2])/5;
   }
	thisbuffer[num_spots-2]=(thisbuffer[num_spots-4]+2*thisbuffer[num_spots-3]+3*thisbuffer[num_spots-2]+4*thisbuffer[num_spots-1])/10;
	thisbuffer[num_spots-1]=(3*thisbuffer[num_spots-1]+2*thisbuffer[num_spots-2]+thisbuffer[num_spots-3]-thisbuffer[num_spots-5])/5;
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
	return(arr1);
}

//smoothing with 3rd order global polynomial fit
/* I'll finish when I have some spare time*/
ArrayStat *Polyfit(ArrayStat *arr1, ArrayStat *arr2){
int i,j,k, row, column, num_spots;
LowessData **ld;
float min, min_diff;
float *buf, *thisbuffer;
int order=3;
float x, y, b[3],a[3][3],c[3],xx[3],h,f,s;

	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;

   for(i=0;i<num_spots;i++){
	   f=1;
   	y=ld[i]->logratio;
      x=ld[i]->intensity;
	   for(j=0;j<2*order-1;j++){
   		if(j>order) {
				c[j]+=f;
	      	f*=x;
	      }else{
		      b[j]+=y;
   		   y*=x;
			}
   	}
   }
   for(i=0;i<order;i++){
   	k=i;
   	for(j=0;j<order;j++){
      	a[i][j]=-a[j][i]/a[i][i];
         for(k=i+1;k<order;k++){
         	a[j][k]+=a[j][i]*a[i][k];
         }
         b[j]+=a[j][i]*b[i];
      }
   }
   xx[order]=b[order]/a[order][order];
   for(i=order-1;i;i--){
   	h=b[i];
      for(j=i+1;j<order;j++){
      	h-=xx[j]*a[i][j];
      }
      xx[i]=h/a[i][i];
   }
   for(i=0;i<num_spots;i++){
   	s=0;
      for(j=order;j<2;j--){
      	s=(s+xx[j])*ld[j]->intensity;
      }
      buf[i]=s+xx[1];
   }

	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
	return(arr1);
}

//Classic LOWESS procedure with parameter - window (percentile)
ArrayStat *Lowess(ArrayStat *arr1, ArrayStat *arr2, int window){
int i,j,k, row, column, flag, num_spots, start, N;
LowessData **ld;
float min, min_diff, prev, a, b, b0, b1, c, d,y;
float *buf, *buf1, *thisbuffer, *thatbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	buf1=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf1[i]=buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thatbuffer=buf1;
   k=0;

   N=(num_spots*window)/100;
   start=0;
   while(start<num_spots-N){
	   a=b=c=d=0;
   	for(i=0;i<N;i++){
         y=thisbuffer[start+i];
		   a=a+i+1; b=b+y; c=c+(i+1)*(i+1); d+=(i+1)*y;
         b1=(a*b-N*d)/(a*a-N*c);
		   b0=(b-b1*a)/N;
      }
   	for(i=0;i<N;i++){
		   thisbuffer[start+i]=b0+b1*i;
      }
      start++;
   }

	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_restorer);
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
   free(buf1);
	return(arr1);
}

ArrayStat *LZscore(ArrayStat *arr1, ArrayStat *arr2, int window){
int i,j,k, row, column, flag, num_spots, start, N;
LowessData **ld;
VectorObject *zz;
float min, min_diff, prev, a, b, b0, b1, c, d, y;
float *buf, *buf1, *thisbuffer, *thatbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	buf1=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf1[i]=buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thatbuffer=buf1;
   k=0;

   N=(int)((num_spots*window)/100.);
   start=0;
   zz=new VectorObject(window,FLOAT);
   a=zz->Median();
	b=zz->StDevMed();
   for(i=0;i<N;i++){
   	zz->vec[i]=buf[i];
   }
   a=zz->Median();
	b=zz->StDevMed();
   for(i=0;i<N/2;i++){
   	zz->vec[i]=(zz->vec[i]-a)/b;
   }
   while(i<num_spots-N/2){
   	for(j=1;j<N;j++) zz->vec[j-1]=zz->vec[j];
      zz->vec[j]=buf[i+N/2];
      buf[i]=(buf[i]-zz->Median())/zz->StDevMed();
      i++;
   }
   while(i<num_spots){
   	zz->vec[i]=(zz->vec[i]-a)/b;
      i++;
   }
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_restorer);
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
   free(buf1);
	return(arr1);
}

/*
ArrayStat *Lowess(ArrayStat *arr1, ArrayStat *arr2, FILE *f){
int i,j,k, row, column, flag, num_spots;
LowessData **ld;
float min, min_diff, prev;
float *buf, *buf1, *thisbuffer, *thatbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	buf1=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf1[i]=buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thatbuffer=buf1;
   k=0;
	   thisbuffer[0]=(39*thisbuffer[0]+8*thisbuffer[1]-4*(thisbuffer[2]+thisbuffer[3]+thisbuffer[4])+thisbuffer[5]-2*thisbuffer[6])/42;
		thisbuffer[1]=(8*thisbuffer[0]+19*thisbuffer[1]+16*thisbuffer[2]+6*thisbuffer[3]-4*thisbuffer[4]-7*thisbuffer[5]+4*thisbuffer[7])/42;
		thisbuffer[2]=(-4*thisbuffer[0]+16*thisbuffer[1]+19*thisbuffer[2]+12*thisbuffer[3]+2*thisbuffer[4]-4*thisbuffer[5]+thisbuffer[6])/42;
	   for(i=3;i<num_spots-3;i++){
	     	thisbuffer[i]=(7*thisbuffer[i]+6*(thisbuffer[i+1]+thisbuffer[i-1])+3*(thisbuffer[i+2]+thisbuffer[i-2])-2*(thisbuffer[i+3]+thisbuffer[i-3]))/21;
	   }
/*
	while(1){
		flag=1;
      k++;
	   for(i=0;i<num_spots-1;i++){
   		if(fabs(thisbuffer[i]-thatbuffer[i])>min_diff){
   		//if(fabs(ld[i]->lr-thisbuffer[i])>min_diff){
         	flag=0;
            break;
			}
      }
      if(flag){
      	break;
      }else{
		   memcpy(thatbuffer,thisbuffer,num_spots*sizeof(float));
         if(k==10)break;
      }
   }

		thisbuffer[num_spots-3]=(thisbuffer[num_spots-7]-4*thisbuffer[num_spots-6]+2*thisbuffer[num_spots-5]+12*thisbuffer[num_spots-4]+19*thisbuffer[num_spots-3]+16*thisbuffer[num_spots-2]-4*thisbuffer[num_spots-1])/42;
		thisbuffer[num_spots-2]=(4*thisbuffer[num_spots-7]-7*thisbuffer[num_spots-6]-4*thisbuffer[num_spots-5]+6*thisbuffer[num_spots-4]+16*thisbuffer[num_spots-3]+19*thisbuffer[num_spots-2]+8*thisbuffer[num_spots-1])/42;
		thisbuffer[num_spots-1]=(-2*thisbuffer[num_spots-7]+4*thisbuffer[num_spots-6]+thisbuffer[num_spots-5]-4*thisbuffer[num_spots-4]-4*thisbuffer[num_spots-3]+8*thisbuffer[num_spots-2]+39*thisbuffer[num_spots-1])/42;
//	middle=logl((arr1->get_trimmed_mean(10)+MINPOS)/(arr2->get_trimmed_mean(10)+MINPOS))/logl(2);
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_restorer);
   lowess_write(ld, num_spots, f);
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
   free(buf1);
	return(arr1);
}
*/



ArrayStat *Slidemean(ArrayStat *arr1, ArrayStat *arr2, FILE *f){
int i,j,k, row, column, flag, num_spots;
LowessData **ld;
float min, min_diff, prev;
float *buf, *buf1, *thisbuffer, *thatbuffer;
	num_spots=arr1->numRows*arr1->numColumns;
   ld=(LowessData **)malloc(num_spots*sizeof(LowessData *));
	buf=(float *)malloc(num_spots*sizeof(float));
	buf1=(float *)malloc(num_spots*sizeof(float));
	k=0;
	for(i=0;i<arr1->numRows;i++){
   	for(j=0;j<arr1->numColumns;j++){
			ld[k]=(LowessData *)malloc(sizeof(LowessData));
      	ld[k]->logratio=ld[k]->lr=logl((arr1->matrix[i][j]+MINPOS)/(arr2->matrix[i][j]+MINPOS))/logl(2);
         ld[k]->intensity=logl(arr1->matrix[i][j]*arr2->matrix[i][j]+MINPOS)/logl(2);
         ld[k]->row=i; ld[k]->column=j; ld[k]->number=k;
         k++;
      }
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_sorter);
	for(i=0;i<num_spots;i++){
   	buf1[i]=buf[i]=ld[i]->logratio;
   }
	min_diff=FLT_MAX;
   for(i=0;i<num_spots-1;i++){
	   for(j=i;i<num_spots-1;i++){
	   	min=fabs(buf[i]-buf[j]);
	      if(min==0) continue;
   	   if(min_diff>min) min_diff=min;
      }
   }
	thisbuffer=buf;
   thatbuffer=buf1;
	while(1){
		thisbuffer[0]=(thisbuffer[0]+thisbuffer[1])/2;
	   for(i=1;i<num_spots-1;i++){
      	thisbuffer[i]=(thisbuffer[i-1]+thisbuffer[i]+thisbuffer[i+1])/3;
      }
      thisbuffer[i]=(thisbuffer[i]+thisbuffer[i-1])/2;
		flag=1;
	   for(i=0;i<num_spots;i++){
   		if(fabs(thisbuffer[i]-thatbuffer[i])>min_diff){
         	flag=0;
            break;
			}
      }
      if(flag){
      	break;
      }else{
		   memcpy(thatbuffer,thisbuffer,num_spots*sizeof(float));
      }
   }
	//middle=logl((arr1->get_trimmed_mean(parameters->trimpercent)+MINPOS)/(arr2->get_trimmed_mean(parameters->trimpercent)+MINPOS))/logl(2);
	for(i=0;i<num_spots;i++){
   	ld[i]->logratio=buf[i];
   }
	for(i=0;i<num_spots;i++){
   	row=ld[i]->row; column=ld[i]->column;
     	arr1->matrix[row][column]=ld[i]->lr - ld[i]->logratio;// + middle;
   }
	qsort(ld, num_spots, sizeof(LowessData **), ld_restorer);
   lowess_write(ld, num_spots, f);
	for(i=0;i<num_spots;i++){
   	free(ld[i]);
   }
   free(ld);
   free(buf);
   free(buf1);
	return(arr1);
}

#define pi 3.141592654

//input vector contains raw measurements
//sinc and cosc are emptry vectors, upon return
//filled with sine and cosin coefficients
//returns vector of periodogram
VectorObject *DirectDFT(VectorObject *input, VectorObject *sinco, VectorObject *cosco){
int i,k,N;
float arg;
VectorObject *pgram;

   N=input->number-input->number%2;
	pgram=new VectorObject(N/2,FLOAT);
	pgram->name=(char *)calloc(strlen(input->name)+16,sizeof(char));
   strcpy(pgram->name,input->name);
   for (i=0;i<N/2;i++) {
      sinco->vec[i] = 0;
      cosco->vec[i] = 0;
      for (k=0;k<N;k++){
         cosco->vec[i] += input->vec[k] * cos((2.0 * k * pi *i)/(float)N);
         sinco->vec[i] += input->vec[k] * sin((2.0 * k * pi *i)/(float)N);
      }
      if(cosco->vec[i]) cosco->vec[i]/=(float)(N/2);
      if(sinco->vec[i]) sinco->vec[i]/=(float)(N/2);
   }

   for(i=0;i<N/2;i++){
   	pgram->vec[i]=sqrt(sinco->vec[i]*sinco->vec[i]+cosco->vec[i]*cosco->vec[i]);
	}
   return(pgram);
}

VectorObject *DirectFFT(VectorObject *input){
//FFT(short int dir,long m,double *x,double *y)
long n,i,i1,j,k,i2,l,l1,l2;
double c1,c2,tx,ty,t1,t2,u1,u2,z, *y, *x;
VectorObject *output;

   /* Calculate the number of points */
   n = 1;
   //for (i=0;i<m;i++) n *= 2;
   n<<input->number;

   output=new VectorObject(input->number,FLOAT);
   memcpy(output->vec,input->vec,input->number*sizeof(float));
   x=(double *)calloc(input->number,sizeof(double));
   y=(double *)calloc(input->number,sizeof(double));

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         //ty = y[i];
         x[i] = x[j];
         //y[i] = y[j];
         x[j] = tx;
         //y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<input->number;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
 //     if (dir == 1)
//         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
//   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
//   }

   memcpy(output->vec,x,input->number*sizeof(float));
   free(y);
   free(x);
   return(output);
}

/*
int m,n,i,i0,i1;
float **x;
VectorObject *spectrum;
	m?
	n=2<<m;
	x=(float *)calloc(n, sizeof(float));
   y=(float *)calloc(n, sizeof(float));
   spectrum=new VectorObject(n);

   for(i=0;i<i0;i++) i+=i0;

   for(l=0;l<m;l++){
   	e=2<<(m+2-l);
      f=e/2;
      u=1;v=0;
      z=3.141838/f;
      c=cos(z);
      s=d*sin(z);
      for(j=0;j<f;j++){
      	for(i=j;i<n;i+=e){
         	o=t*u+r*v;
            x[i-1]=p;
            y[i-1]=q;
         }
         w==u*c-v*s;
         v=v*c+u*s;
         u=w;
      }
   }
   j=1;
   for(i=0;i<n-1;i++){
   	if(i>=j) goto:l150;
      	j1=j-1;
         i1=i-1;
         p=x[j1];
         q=y[j1];
         x[j1]=x[i1];
         y[j1]=y[i1];
         x[i1]=p;
         y[i1]=q;

      l150:k=n/2;
      l160:if(k>=j)goto l180;
      	j=j-k;
         k=k/2;
         goto l160;
      l180:j=j+k;
   }
   for(k=0;k<n-1;k++){
   	a=sqrt(x[k]*x[k]+y[k]*y[k]);
   	q=0;
      if(a==0) goto l270;
      q=acos(x[k]/a);
      if(y[k]<0 then q=-q;



   free(x);
   free(y);
   return(spectrum);
} */
