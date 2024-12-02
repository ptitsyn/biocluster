//Copyright(c) Andrey Ptitsyn, PBRC 2002
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>


#define CH1INTENSITY 1
#define CH2INTENSITY 2
#define CH1BACKGROUND 3
#define CH2BACKGROUND 4
#define CH1INTENSITYSTDEV 5
#define CH2INTENSITYSTDEV 6
#define CH1BACKGROUNDSTDEV 7
#define CH2BACKGROUNDSTDEV 8
#define CH1DIAMETER 9
#define CH2DIAMETER 10
#define CH1AREA 11
#define CH2AREA 12
#define CH1FOOTPRINT 13
#define CH2FOOTPRINT 14
#define CH1CIRCULARITY 15
#define CH2CIRCULARITY 16
#define CH1SPOTUNIFORMITY 17
#define CH2SPOTUNIFORMITY 18
#define CH1BKGUNIFORMITY 19
#define CH2BKGUNIFORMITY 20
#define CH1SIGNALNOISERATIO 21
#define CH2SIGNALNOISERATIO 22
#define CH1CONFIDENCE 23
#define CH2CONFIDENCE 24

#define ATLAS7909_1 0
#define PLASTIC_RAT4K 0

#define GENCODE 0
#define INTENSITY_1 1
#define BACKGROUND_1 2
#define ADJINTENSITY_1 3
#define INTENSITY_2 4
#define BACKGROUND_2 5
#define ADJINTENSITY_2 6
#define RATIO 7
#define DIFFERENCE 8
#define PROTEINGENE 9

#define GENERIC 0
#define SCANARRAY 1
#define GENEPIX 2


#define MINPOS 0.000001

//////////////////////////////////////////////////////////////////////////////

class ScanArraySpot{
public:
	int Number;
   int ArrayRow;
   int ArrayColumn;
   int Row;
   int Column;
   char *Name;
   int XLocation;
   int YLocation;
   float ch1Intensity;
   float ch1Background;
   float ch1IntensityStdDev;
   float ch1BackgroundStdDev;
   float ch1Diameter;
   int ch1Area;
   float ch1Footprint;
   float ch1Circularity;
   float ch1SpotUniformity;
   float ch1BkgUniformity;
   float ch1SignalNoiseRatio;
   float ch1Confidence;
   float ch2Intensity;
   float ch2Background;
   float ch2IntensityStdDev;
   float ch2BackgroundStdDev;
   float ch2Diameter;
   int ch2Area;
   float ch2Footprint;
   float ch2Circularity;
   float ch2SpotUniformity;
   float ch2BkgUniformity;
   float ch2SignalNoiseRatio;
   float ch2Confidence;
   short int ch1DiameterFilter;
   short int ch1AreaFilter;
   short int ch1FootprintFilter;
   short int ch1CircularityFilter;
   short int ch1SpotUniformityFilter;
   short int ch1BkgUniformityFilter;
   short int ch1SignalNoiseRatioFilter;
   short int ch1ReplicateFilter;
   short int ch2DiameterFilter;
   short int ch2AreaFilter;
   short int ch2FootprintFilter;
   short int ch2CircularityFilter;
   short int ch2SpotUniformityFilter;
   short int ch2BkgUniformityFilter;
   short int ch2SignalNoiseRatioFilter;
   short int ch2ReplicateFilter;
   short int IgnoreFilter;

   ScanArraySpot();
   ScanArraySpot(FILE *, int num_data_columns);
   ~ScanArraySpot();

   int Read(FILE *f, int num_data_columns);
   int Write(FILE *f);
   int Write(FILE *f, float extravalue);
   int Write(FILE *f, float extravalue1, float extravalue2, float extravalue3);
   int SummarizeFilters();
   void *ValueOf(int code);
   void FillFrom(ScanArraySpot *);
};

class ScanArrayHeader{
public:
	char *UserName;
	char *Computer;
	char *Date;
	char *Experiment;
	char *ExperimentPath;
	char *Protocol;
	int Version;
   ScanArrayHeader();
   ~ScanArrayHeader();
	int Find(FILE *);
   int Read(FILE *);
   int Write(FILE *f);
};

class ProtocolInfo{
public:
	char *Units;
	int ArrayRows;
	int ArrayColumns;
	int Rows;
	int Columns;
	float ArrayRowSpacing;
	float ArrayColumnSpacing;
	int SpotRowSpacing;
	int SpotColumnSpacing;
	int SpotDiameter;
	int Interstitial;
   char *InterstitialComment;
	int SpotsPerArray;
	int TotalSpots;
	int CrosstalkCorrected;
	int BackgroundSubtracted;
	char *QuantificationMethod;
	char *QualityConfidenceCalculation;
   ProtocolInfo();
   ~ProtocolInfo();
   int Read(FILE *);
   int Write(FILE *f);
   int Find(FILE *);
};

class ToleranceWeight{
public:
	float DiameterMin;
	float DiameterMax;
	float DiameterWeight;
	float AreaMin;
	float AreaMax;
	float AreaWeight;
	float FootprintMin;
	float FootprintMax;
	float FootprintWeight;
	float CircularityMin;
	float CircularityMax;
	float CircularityWeight;
	float SpotUniformityMin;
	float SpotUniformityMax;
	float SpotUniformityWeight;
	float BkgUniformityMin;
	float BkgUniformityMax;
	float BkgUniformityWeight;
	float SignalNoiseRatioMin;
	float SignalNoiseRatioMax;
	float SignalNoiseRatioWeight;
	float ReplicateUniformityMin;
	float ReplicateUniformityMax;
	float ReplicateUniformityWeight;

   ToleranceWeight();
   ~ToleranceWeight();
	int Find(FILE *);
   int Read(FILE *);
   int Write(FILE *f);
};

typedef struct{
	char *ImageFile;
   char *Fluorophor;
   char *Barcode;
   char *Units;
   int	XUnitsPerPixel;
   int	YUnitsPerPixel;
   int	XOffset;
   int	YOffset;
   char *Status;
}CHANNEL;

class ImageInfo{
public:
	CHANNEL ch1;
	CHANNEL ch2;
	ImageInfo();
   ~ImageInfo();
	int Find(FILE *);
	int Read(FILE *);
   int Write(FILE *f);
};


class ScanArraySlideData{
public:
	int num_spots;
   int num_arrays;
   int num_subarrays;
   int num_array_columns;
   int num_array_raws;
   int num_data_columns;
   ScanArraySpot **spot;
   ScanArrayHeader *header;
   ProtocolInfo *protocol;
   ToleranceWeight *tolerance;
   ImageInfo *image;

   ScanArraySlideData();
   ~ScanArraySlideData();
   int Read(FILE *);
   int Find(FILE *);
   int Write(FILE *f);
   int WriteShort(FILE *f);
   int WriteShort(FILE *f, char *extraname, float *extravalue);
   int WriteShort(FILE *f, char *extraname1, char *extraname2, char *extraname3, float *extravalue1, float *extravalue2, float *extravalue3);
   int ReturnValuesFrom(class ArrayStat *stat, int code);
   int SummarizeFilters();
};

class ArrayStat:public MatrixObject{
public:
	int code;
	int subarray;
	int offset;
	float Mean;
   float Median;
   float Stdev;
   float MAD;
   float TMAD;
   ArrayStat();
   ArrayStat(ArrayStat *);
   ArrayStat(ScanArraySlideData *);
   ArrayStat(ScanArraySlideData *, int code);
   ArrayStat(ScanArraySlideData *, int field, int arraynumber);
   ArrayStat(ScanArraySlideData *, int code, int subarray_col, int subarray_raw);
   ~ArrayStat();

   float get_mean();
   float get_median();
   float get_MAD();
   float get_TMAD(int percentile);
   float get_trimmed_mean(int percentile);
   float get_stdev();
   void ScaleToMean();
   void ScaleToMedian();
   void ScaleToMAD();
   void ScaleToTMAD(int percentile);
   void ScaleToStdev();
   void ScaleTo(float factor);
   void ZScore();
   void GetLog2();
   void shuffle_by_xy(int times);
   int mean_equal_to(ArrayStat *other, int bon, float *d);
   int grid_is_uniform( int grid );
   int grid_is_bonferroni_uniform( int grid );
	int Write(FILE *);

};

class ArrayFile{
public:
	char *name;
   FILE *f;
};

typedef struct{
	float logratio;
   float lr;
   float intensity;
   int row;
   int column;
   int number;
}LowessData;

VectorObject *GetSlideControls(int code, char *keyword, ScanArraySlideData *slide);

VectorObject *GetSlideSubarrayControls(int code, char *keyword, ScanArraySlideData *slide, int ArrayRow, int ArrayColumn);

int lowess_write(LowessData **ld, int num_spots, FILE *f);

ArrayStat *Slidemean(ArrayStat *arr1, ArrayStat *arr2);

ArrayStat *Slidemean(ArrayStat *arr1, ArrayStat *arr2, FILE *f);

ArrayStat *Polysmooth(ArrayStat *arr1, ArrayStat *arr2);

ArrayStat *Lsmooth(ArrayStat *arr1, ArrayStat *arr2);

ArrayStat *Polyfit(ArrayStat *arr1, ArrayStat *arr2);

ArrayStat *Lowess(ArrayStat *arr1, ArrayStat *arr2, int window);

ArrayStat *LZscore(ArrayStat *arr1, ArrayStat *arr2, int window);

VectorObject *DirectDFT(VectorObject *input, VectorObject *sinco, VectorObject *cosco);

VectorObject *DirectFFT(VectorObject *input);

////////////////////////////////////////////////////////////////////////////////

class AtlasArraySpot{
public:
	char *Genecode;
	int Intensity_1;
	int Background_1;
	int AdjIntensity_1;
	int Intensity_2;
	int Background_2;
	int AdjIntensity_2;
	float Ratio;
	int Difference;
	char *ProteinGene;
   int IgnoreFilter;

   AtlasArraySpot();
   AtlasArraySpot(FILE *);
   ~AtlasArraySpot();

   int Read(FILE *f);
   int Write(FILE *f);
   int Write(FILE *f, float extravalue);
   void *ValueOf(int code);
   void FillFrom(AtlasArraySpot *);
};

class AtlasArraySlideData{
public:
	int num_spots;
   int num_arrays;
   int num_subarrays;
   int num_array_columns;
   int num_array_rows;
   int num_data_rows;
   int num_data_columns;
   AtlasArraySpot **spot;
   char *header;

   AtlasArraySlideData();
   ~AtlasArraySlideData();
   int Read(FILE *);
   int Read(FILE *, int layout);
   int Find(FILE *);
   int Write(FILE *f);
   int Write(FILE *f, char *extraname, float *extravalue);
   int ReturnValuesFrom(class ArrayStat *stat, int code);
   int SummarizeFilters();
   ArrayStat *MakeStat(int layout, int code);
};

////////////////////////////////////////////////////////////////////

class GenePixArraySpot{
};

class GenePixArraySlideData{
};
////////////////////////////////////////////////////////////////////

class GenericArraySpot{
public:
	int Number;
   int ArrayRow;
   int ArrayColumn;
   int Row;
   int Column;
	char *Name;
	int Intensity_1;
	int Background_1;
	int Intensity_2;
	int Background_2;
	float Ratio;
	int Difference;
	char *ProteinGene;
   int IgnoreFilter;

   GenericArraySpot();
   GenericArraySpot(FILE *);
   ~GenericArraySpot();

   int Read(FILE *f);
   int Write(FILE *f);
   int Write(FILE *f, float extravalue);
   void *ValueOf(int code);
   void FillFrom(AtlasArraySpot *);
};

class GenericArraySlideData{
public:
	int num_spots;
   int num_arrays;
   int num_subarrays;
   int num_array_columns;
   int num_array_rows;
   int num_data_rows;
   int num_data_columns;
   int ch1intensity;
   int ch1background;
   int ch2intensity;
   int ch2background;
   GenericArraySpot **spot;
   char *header;

   GenericArraySlideData();
   ~GenericArraySlideData();
   int Read(FILE *);
   int Read(FILE *, int layout);
   int Find(FILE *);
   int Write(FILE *f);
   int Write(FILE *f, char *extraname, float *extravalue);
   int ReturnValuesFrom(class ArrayStat *stat, int code);
	int ReadLayout(FILE *f);
   ArrayStat *MakeStat(int field);
   ArrayStat *MakeStat(int field, int arraynumber);
   ArrayStat *MakeStat(int field, int subarray_row, int subarray_col);
};


