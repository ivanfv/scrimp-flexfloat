/* #############################################################################

	         ███████╗ ██████╗██████╗ ██╗███╗   ███╗██████
                 ██╔════╝██╔════╝██╔══██╗██║████╗ ████║██╔══██╗
                 ███████╗██║     ██████╔╝██║██╔████╔██║██████╔╝
                 ╚════██║██║     ██╔══██╗██║██║╚██╔╝██║██╔═══╝ 
                 ███████║╚██████╗██║  ██║██║██║ ╚═╝ ██║██║ 
	         ╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     
                 +-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+
                 | F l e x F l o a t        B e n c h m a r k |
                 +-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+

SCRIMP FlexFloat benchmark provides a precision exploation tool for time series
analysis using Matrix Profile. This program allows changing the number of bits
used by floating point operations to meassure how affects to the accuracy of 
result. This program is not inteded to be used for performance evaluation
purposes as long as it is based on a floating point emulation library.

The code is based on the SCRIMP C++ implementation developed by Yan Zhu, 
Chin-Chia Michael Yeh, Zachary Zimmerman, Kaveh Kamgar and Eamonn Keogh that 
can be found at https://sites.google.com/site/scrimpplusplus/. It was modified
by Ivan Fernandez in August 2019 to include the FlexFloat library.

Details of the SCRIMP++ algorithm can be found at:
    - "SCRIMP++: Motif Discovery at Interactive Speeds", ICDM 2018.
    - Matrix Profile website: https://www.cs.ucr.edu/~eamonn/MatrixProfile.html

Details of Flexfloat library can be found at:
    - "FlexFloat: A Software Library for Transprecision Computing",
             IEEE Transactions on Computer-Aided Design of Int. Cir. and Systems
    - FlexFloat repository: https://github.com/oprecomp/flexfloat    

********************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
********************************************************************************
Usage: 
>> scrimpplusplus InputFile WindowSize nThreads Scale
	- InputFile: Name of the time series file
	- WindowSize: Subsequence length m
	- nThreads: Number of threads to be spawn
	- Scale: Scale factor for the time series data

Example input:
>> scrimp_ff taxi.txt 100 8 1

Example output:
The code will generate one output file (e.g. result_taxi.csv)
	- First column is the index of the data.
	- Second column is the matrix profile flexfloat distance value.
	- Third column is the matrix profile flexfloat index.
	- Fourth column is the matrix profile distance value.
	- Fifth column is the matrix profile index.
	- Sixth column is the error
############################################################################# */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <omp.h>
#include "../flexfloat/include/flexfloat.h"

#define PATH_TSERIES "./timeseries/"
#define PATH_CFG "./configs/"
#define PATH_RESULT "./results/result_"
#define EXCLUSION_FACTOR 4

unsigned dist_exp;
unsigned dist_man;
unsigned prof_exp;
unsigned prof_man;
unsigned dotp_exp;
unsigned dotp_man;;
unsigned stats_exp;
unsigned stats_man;

static struct timeval tm1;

void scrimp(double * tSeries,  double * AMean, double * ASigma, 
		int timeSeriesLength, int ProfileLength,
		int windowSize, int* idx, double * profile, int * profileIndex,
	       	int exclusionZone, int numThreads)
{
	/* Private structures initialization -------------------------------- */
	double * profile_tmp   = malloc(sizeof(double) 
			* ProfileLength * numThreads);
	int * profileIndex_tmp = malloc(sizeof(int)    
			* ProfileLength * numThreads);
	for(int i = 0; i < ProfileLength * numThreads; i++)
	{
		profile_tmp[i]      = INFINITY;
		profileIndex_tmp[i] = 0;
	}
	/* ------------------------------------------------------------------ */

	#pragma omp parallel
	{
		double distance, windowSizeDTYPE, lastz;

		int diag, my_offset, i, j, ri;

		windowSizeDTYPE = (double) windowSize;

		my_offset = omp_get_thread_num() * ProfileLength;

		#pragma omp for schedule(dynamic)
		for (ri = 0; ri < (ProfileLength - (exclusionZone + 1)); ri++)
		{
			diag = idx[ri];
			lastz = 0;

			/* Dot product calculation -------------------------- */
			for (j = diag; j < windowSize + diag; j++)
			{
				lastz += tSeries[j] * tSeries[j - diag];
			}

			j = diag;
			i = 0;

			/* Distance calculation ----------------------------- */
			distance = 2 * (windowSizeDTYPE - (lastz - 
				windowSizeDTYPE * AMean[j] * AMean[i]) / 
				(ASigma[j] * ASigma[i]));

			/* -------------------------------------------------- */
			
			/* Profile update ----------------------------------- */
			if (distance < profile_tmp[my_offset + j])
			{
				profile_tmp[my_offset + j]     = distance;
				profileIndex_tmp [my_offset+j] = i;
			}
			if (distance < profile_tmp[my_offset + i])
			{
				profile_tmp[my_offset + i]       = distance;
				profileIndex_tmp [my_offset + i] = j;
			}
			/* -------------------------------------------------- */
			i = 1;

			for(j = diag + 1; j< ProfileLength; j++)
			{
				/* Dot product update ----------------------- */
				lastz += (tSeries[j + windowSize - 1] * 
					tSeries[i + windowSize - 1]) - 
					(tSeries[j - 1] * tSeries[i - 1]);
				/* ------------------------------------------ */

				/* Distance calculation --------------------- */
				distance =  2 * (windowSizeDTYPE - (lastz -  
					AMean[j]  * AMean[i] * windowSizeDTYPE) 
					/ (ASigma[j] * ASigma[i]));
				/* ------------------------------------------ */

				/* Profile update --------------------------- */
				if (distance < profile_tmp[my_offset + j])
				{
					profile_tmp[my_offset + j] = distance;
					profileIndex_tmp [my_offset+ j] = i;
				}

				if (distance < profile_tmp[my_offset + i])
				{
					profile_tmp[my_offset + i] = distance;
					profileIndex_tmp[my_offset + i] = j;
				}
				/* ------------------------------------------ */
				i++;
			}

		}

		#pragma omp barrier
		
		/* Final profile reduction ---------------------------------- */ 
		double min_distance;
		int min_index;

		#pragma omp for schedule(static)
		for (int colum = 0; colum < ProfileLength; colum++)
		{
			min_distance = INFINITY;

			min_index = 0;

			for(int row = 0; row < numThreads; row++)
			{
				if(profile_tmp[colum + (row*ProfileLength)] 
						< min_distance)
				{
					min_distance = profile_tmp[colum + 
						(row * ProfileLength)];
					min_index    = profileIndex_tmp[colum 
						+ (row * ProfileLength)];
				}
			}
			profile[colum]      = min_distance;
			profileIndex[colum] = min_index;
		}
		#pragma omp barrier
		/* ---------------------------------------------------------- */
	}
	free(profile_tmp);
	free(profileIndex_tmp);
}

void scrimp_ff(flexfloat_t * tSeries, flexfloat_t* AMean, flexfloat_t* ASigma, 
		int timeSeriesLength,int ProfileLength,	flexfloat_t windowSize, 
		int* idx, flexfloat_t * profile, int* profileIdxs, 
		int exclusionZone, int numThreads)
{

	int win = (int) ff_get_double(&windowSize);

	/* Private structures initilization --------------------------------- */
	flexfloat_t * profile_priv;
	int         * profileIdxs_priv;
	profile_priv    = malloc(sizeof(flexfloat_t) 
			* timeSeriesLength * numThreads);
	profileIdxs_priv = malloc(sizeof(int) 
			* timeSeriesLength * numThreads);
	for(int i = 0; i < timeSeriesLength * numThreads; i++)
	{
		ff_init_double(&profile_priv[i], INFINITY, 
				(flexfloat_desc_t) {prof_exp, prof_man});
		profileIdxs_priv[i] = 0;
	}
	/* ------------------------------------------------------------------ */

	#pragma omp parallel
	{
		flexfloat_t substr;
		flexfloat_t distance;
		flexfloat_t sigma_prods;
		flexfloat_t mean_prods;
		flexfloat_t constant_2;
		flexfloat_t lastz;
		flexfloat_t lastz_cast;
		flexfloat_t dist_cast;
		flexfloat_t mean_cast;
		flexfloat_t sigma_cast;
		int myoffset = omp_get_thread_num() * timeSeriesLength;

		ff_init_double(&substr, 0, 
				(flexfloat_desc_t) {dotp_exp, dotp_man});
		ff_init_double(&distance, 0, 
				(flexfloat_desc_t) {dist_exp, dist_man});
		ff_init_double(&mean_prods, 0, 
				(flexfloat_desc_t) {stats_exp, stats_man});
		ff_init_double(&sigma_prods, 0, 
				(flexfloat_desc_t) {stats_exp, stats_man});
		ff_init_double(&constant_2, 2.0, 
				(flexfloat_desc_t) {dist_exp, dist_man});

		#pragma omp for schedule(dynamic)
		for (int ri = 0; ri < (ProfileLength - 
					(exclusionZone + 1)); ri++)
		{
			int subseq = idx[ri];

			/* Dot product calculation -------------------------- */
			ff_init_double(&lastz,0, (flexfloat_desc_t) 
					{dotp_exp, dotp_man});
			for (int w = 0; w < win; w++)
			{
				ff_fma(&lastz,  &tSeries[w + subseq], 
						&tSeries[w], &lastz);
			}
			ff_cast(&lastz_cast, &lastz, (flexfloat_desc_t) 
					{dist_exp, dist_man});
			/* -------------------------------------------------- */

			/* Distance calculation ----------------------------- */
                        ff_mul(&sigma_prods, &ASigma[subseq], &ASigma[0]);
			ff_mul(&mean_prods, &AMean[subseq], &AMean [0]);
			ff_cast(&mean_cast, &mean_prods, (flexfloat_desc_t)
                                        {dist_exp, dist_man});
                        ff_cast(&sigma_cast, &sigma_prods, (flexfloat_desc_t)
					{dist_exp, dist_man});
   			ff_mul(&distance,    &mean_cast,   &windowSize);
                        ff_sub(&distance,    &lastz_cast, &distance);
                        ff_div(&distance,    &distance,   &sigma_cast);
                        ff_sub(&distance,    &windowSize, &distance);
                        ff_mul(&distance,    &distance,   &constant_2);
			/* -------------------------------------------------- */

			/* Profile update ----------------------------------- */
			ff_cast(&dist_cast, &distance, (flexfloat_desc_t)
                                        {prof_exp, prof_man});
			if (ff_lt(&dist_cast, &profile_priv[subseq + myoffset]))
			{
				profile_priv[subseq + myoffset]     = dist_cast;
				profileIdxs_priv[subseq + myoffset] = 0;
			}
			if (ff_lt(&dist_cast, &profile_priv[myoffset]))
			{
				profile_priv[myoffset]     = dist_cast;
				profileIdxs_priv[myoffset] = subseq;
			}
			/* -------------------------------------------------- */      
			int i = 1;

			for(int j = subseq + 1; j< ProfileLength; j++)
			{
				/* Dot product update ----------------------- */
				ff_fma(&lastz,  &tSeries[j + win - 1], 
						&tSeries[i + win - 1],  &lastz);
				ff_mul(&substr, &tSeries[j - 1],        
						&tSeries[ i - 1]);
				ff_sub(&lastz,  &lastz, &substr);
				ff_cast(&lastz_cast, &lastz, (flexfloat_desc_t)
					       	{dist_exp, dist_man});
				/* ------------------------------------------ */

				/* Distance calculation --------------------- */
				ff_mul(&sigma_prods, &ASigma[j],  &ASigma[i]);
				ff_mul(&mean_prods,  &AMean[j],   &AMean [i]);
                        	ff_cast(&mean_cast, &mean_prods, 
					(flexfloat_desc_t){dist_exp, dist_man});
                        	ff_cast(&sigma_cast, &sigma_prods, 
					(flexfloat_desc_t){dist_exp, dist_man});
				ff_mul(&distance,    &mean_cast,  &windowSize);
				ff_sub(&distance,    &lastz_cast, &distance);
				ff_div(&distance,    &distance,   &sigma_cast);
				ff_sub(&distance,    &windowSize, &distance);
				ff_mul(&distance,    &distance,   &constant_2);
				/* ------------------------------------------ */

				/* Profile update --------------------------- */
				ff_cast(&dist_cast, &distance,(flexfloat_desc_t)
                                        {prof_exp, prof_man});
				if (ff_lt(&dist_cast,&profile_priv[j+myoffset]))
				{
					profile_priv[j + myoffset]  = dist_cast;
					profileIdxs_priv[j+myoffset] = i;
				}
				if (ff_lt(&dist_cast,&profile_priv[i+myoffset]))
				{
					profile_priv[i + myoffset]  = dist_cast;
					profileIdxs_priv[i+myoffset] = j;
				}
				/* ------------------------------------------ */
				i++;
			}

		}
		#pragma omp barrier

		/* Final profile reduction ---------------------------------- */
		flexfloat_t min_distance;
		#pragma omp for schedule(static)
		for (int colum = 0; colum < ProfileLength; colum++)
		{
			ff_init_double(&min_distance,  INFINITY, 
				(flexfloat_desc_t) {prof_exp,prof_man});
			int min_index = 0;

			for(int row = 0; row < numThreads; row++)
			{
				if(ff_lt(&profile_priv[colum + 
					(row*timeSeriesLength)], &min_distance))
				{
					min_distance = profile_priv[colum + 
						(row *	timeSeriesLength)];
					min_index    = profileIdxs_priv[colum +
					       	(row * timeSeriesLength)];
				}
			}
			profile[colum]     = min_distance;
			profileIdxs[colum] = min_index;
		}
		#pragma omp barrier
		/* ---------------------------------------------------------- */
	}
	free(profile_priv);
	free(profileIdxs_priv);
}

static inline void start()
{
	/* Starting chronograph --------------------------------------------- */
	gettimeofday(&tm1, NULL);
	/* ------------------------------------------------------------------ */
}

static inline void stop()
{
	/* Stop chronograph and display time -------------------------------- */
	struct timeval tm2;
	gettimeofday(&tm2, NULL);
	unsigned long long t = 1000 * (tm2.tv_sec - tm1.tv_sec) + 
		(tm2.tv_usec - tm1.tv_usec) / 1000;
	double seconds = (double) t / 1000;
	printf("[INFO] DONE (elapsed time %.3f seconds)\n", seconds);
	/* ------------------------------------------------------------------ */
}

void read_config(char * file_name)
{
	
	char * path_config = malloc(1000 * sizeof(char));
	
        strcpy(path_config, PATH_CFG);
        strcat(path_config, file_name);
        path_config[strlen(path_config)-4] = 0;
        strcat(path_config, ".cfg");

        FILE * file = fopen (path_config, "r");
	int err1, err2, err3, err4, err5, err6, err7, err8;
	
	if (file == NULL) printf("CFG FILE ERRROR\n");
	err1 = fscanf(file, "%d", &dist_exp);
	err2 = fscanf(file, "%d", &dist_man);
  	err3 = fscanf(file, "%d", &dotp_exp);
   	err4 = fscanf(file, "%d", &dotp_man);
   	err5 = fscanf(file, "%d", &stats_exp);
        err6 = fscanf(file, "%d", &stats_man);
        err7 = fscanf(file, "%d", &prof_exp);
        err8 = fscanf(file, "%d", &prof_man);
	
       if (!err1 || !err2 || !err3 || !err4 || !err5
		      || !err6 || !err7 || !err8) 
	       printf("CFG FILE ERRROR\n");

	free(path_config);
	fclose (file);
}



void print_header()
{
	/* Printing program header ------------------------------------------ */
	printf("\n");	
	printf("███████╗ ██████╗██████╗ ██╗███╗   ███╗██████╗ \n"
	       "██╔════╝██╔════╝██╔══██╗██║████╗ ████║██╔══██╗\n"
	       "███████╗██║     ██████╔╝██║██╔████╔██║██████╔╝\n"
	       "╚════██║██║     ██╔══██╗██║██║╚██╔╝██║██╔═══╝ \n"
	       "███████║╚██████╗██║  ██║██║██║ ╚═╝ ██║██║     \n"
	       "╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     \n"
	       "+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+\n"
	       "| F l e x F l o a t        B e n c h m a r k |\n"    
	       "+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+\n");
	printf("\n");
	/* ------------------------------------------------------------------ */
}

int main(int argc, char* argv[])
{
	FILE   * fp;
	time_t t;
	int    windowSize, ProfileLength, exclusionZone;
	int    timeSeriesLength, numThreads;
	double tSeriesMin, tSeriesMax;
	double elapsed_time, scaleFactor;
	double      * tSeries;
	double      * AMean;
	double      * ASigma;
	double      * profile;
	int         * profileIdxs;
	int         * profileIdxs_ff;
	int         * idx;
	char        * path_tSeries;
	char	    * path_result;
	flexfloat_t * tSeries_ff;
	flexfloat_t * AMean_ff;
	flexfloat_t * ASigma_ff;
	flexfloat_t * profile_ff;

        ff_start_stats();
	print_header();

	/* Getting program arguments ---------------------------------------- */
	if(argc !=5)
	{
		printf("[ERROR] usage: ./scrimp timeseries.txt window_size"
				" num_threads scale_factor\n");
		return 0;
	}
	windowSize  = atoi(argv[2]);
	numThreads  = atoi(argv[3]);
	scaleFactor = atof(argv[4]);
	omp_set_num_threads(numThreads);
	read_config(argv[1]);
	/* ------------------------------------------------------------------ */

	/* Time series loading and memory allocating ------------------------ */	
	path_tSeries = malloc(1000 * sizeof(char));
	strcpy(path_tSeries, PATH_TSERIES);
	strcat(path_tSeries, argv[1]);
	printf("##############################################\n");
	printf("[INFO] Loading %s ...\n", argv[1]);
	timeSeriesLength = 0; 
	fp = fopen(path_tSeries, "r");
	char c;
	for (c = getc(fp); c != EOF; c = getc(fp)) 
		if (c == '\n') 
			timeSeriesLength = timeSeriesLength + 1;
	fclose(fp);

	ProfileLength = timeSeriesLength - windowSize + 1;
	exclusionZone = windowSize / EXCLUSION_FACTOR;

	tSeries          = malloc(sizeof(double) * timeSeriesLength);
	AMean            = malloc(sizeof(double) * ProfileLength);
	ASigma           = malloc(sizeof(double) * ProfileLength);
	profile          = malloc(sizeof(double) * ProfileLength);
	profileIdxs      = malloc(sizeof(int)    * timeSeriesLength);
	profileIdxs_ff   = malloc(sizeof(int)    * timeSeriesLength);
	idx              = malloc(sizeof(int)    * timeSeriesLength);
	tSeries_ff       = malloc(sizeof(flexfloat_t) * timeSeriesLength);
	AMean_ff         = malloc(sizeof(flexfloat_t) * timeSeriesLength);
	ASigma_ff        = malloc(sizeof(flexfloat_t) * timeSeriesLength);
	profile_ff       = malloc(sizeof(flexfloat_t) * timeSeriesLength);

	tSeriesMin = INFINITY;
	tSeriesMax = 0;

	int err;
	fp = fopen(path_tSeries, "r");
	for(int i = 0; i < timeSeriesLength; i++)
	{
		err = fscanf(fp,"%lf",&tSeries[i]);
		if(err < 0) return -1;
		tSeries[i] *= scaleFactor;
		if(tSeries[i] < tSeriesMin) tSeriesMin = tSeries[i];
		if(tSeries[i] > tSeriesMax) tSeriesMax = tSeries[i];
		ff_init_double(&tSeries_ff[i], tSeries[i],
                                (flexfloat_desc_t) {dotp_exp, dotp_man});
	}

	fclose(fp);
	printf("[INFO] DONE\n");
	/* ------------------------------------------------------------------ */

	printf("[INFO] Program parameters:\n");
	printf("----------------------------------------------\n");
	printf("  Time series length: %d\n", timeSeriesLength);
	printf("  Window size:        %d\n", windowSize);
	printf("  Time series min:    %f\n", tSeriesMin);
	printf("  Time series max:    %f\n", tSeriesMax);
	printf("  Dot product max:    %f\n", pow(tSeriesMax,2) * windowSize);
	printf("  Number of threads:  %d\n", numThreads);
        printf("  Scale factor:       %.4f\n", scaleFactor);
        printf("  FF dist - exp, man: %d, %d\n", dist_exp, dist_man);
	printf("  FF dotp - exp, man: %d, %d\n", dotp_exp, dotp_man);
        printf("  FF stat - exp, man: %d, %d\n", stats_exp, stats_man);
        printf("  FF prof - exp, man: %d, %d\n", prof_exp, prof_man);
	printf("----------------------------------------------\n");

	/* Preprocessing statistics ----------------------------------------- */
	printf("[INFO] Preprocessing statistics ...\n");
	double * ACumSum   = malloc(sizeof(double) * timeSeriesLength);
	double * ASum      = malloc(sizeof(double) * timeSeriesLength);
	double * ASumSq    = malloc(sizeof(double) * timeSeriesLength);
	double * ASigmaSq  = malloc(sizeof(double) * timeSeriesLength);
	double * ASqCumSum = malloc(sizeof(double) * timeSeriesLength);

	ACumSum[0] = tSeries[0];
	for (int i = 1; i < timeSeriesLength; i++)
		ACumSum[i] = tSeries[i] + ACumSum[i - 1];  
	ASqCumSum[0] = tSeries[0] * tSeries[0];
	for (int i = 1; i < timeSeriesLength; i++)
		ASqCumSum[i] = tSeries[i] * tSeries[i] + ASqCumSum[i - 1];
	ASum[0] = ACumSum[windowSize - 1];
	for (int i = 0; i < timeSeriesLength - windowSize; i++)
		ASum[i + 1] = ACumSum[windowSize + i] - ACumSum[i];
	ASumSq[0] = ASqCumSum[windowSize - 1];
	for (int i = 0; i < timeSeriesLength - windowSize; i++)
		ASumSq[i + 1] = ASqCumSum[windowSize + i] - ASqCumSum[i];
	for (int i = 0; i < ProfileLength; i++)
		AMean[i] = ASum[i] / windowSize;
	for (int i = 0; i < ProfileLength; i++)
		ASigmaSq[i] = ASumSq[i] / windowSize - AMean[i] * AMean[i];
	for (int i = 0; i < ProfileLength; i++)
		ASigma[i] = sqrt(ASigmaSq[i]);

	free(ACumSum);
	free(ASum);
	free(ASqCumSum);
	free(ASumSq);
	free(ASigmaSq);

	printf("[INFO] DONE\n");
	/* ------------------------------------------------------------------ */

	/* Initializing diagonals and flexfloat statistics------------------- */
	for (int i = exclusionZone + 1; i < ProfileLength; i++)
	{
		idx[i - (exclusionZone + 1)] = i;
	}

	for (int i = 0; i < ProfileLength; i++)
	{
		profile[i]     = INFINITY;
		profileIdxs[i] = 0;
	}

	for(int i = 0; i < ProfileLength; i++)
	{
		ff_init_double(&AMean_ff[i],   AMean[i],   (flexfloat_desc_t) 
				{stats_exp, stats_man});
		ff_init_double(&ASigma_ff[i],  ASigma[i],  (flexfloat_desc_t) 
				{stats_exp, stats_man});
		ff_init_double(&profile_ff[i], profile[i], (flexfloat_desc_t) 
				{prof_exp, prof_man});
	}

	flexfloat_t  windowSize_ff;
	ff_init_double(&windowSize_ff, windowSize, (flexfloat_desc_t) 
			{dist_exp, dist_man});
	/* ------------------------------------------------------------------ */

	/* Running SCRIMP FF ------------------------------------------------ */
	printf("[INFO] Running SCRIMP FlexFloat ...\n");
	start();
	scrimp_ff(tSeries_ff, AMean_ff, ASigma_ff, timeSeriesLength,
			ProfileLength, windowSize_ff, idx,
			profile_ff, profileIdxs_ff, exclusionZone, numThreads);
	stop();
	/* ------------------------------------------------------------------ */

	/* Running SCRIMP double precision ---------------------------------- */
	printf("[INFO] Running SCRIMP  ...\n");

	start();
	scrimp(tSeries, AMean, ASigma, timeSeriesLength, ProfileLength,
			windowSize, idx, profile, profileIdxs, exclusionZone,
			numThreads);
	stop();
	/* ------------------------------------------------------------------ */

	/* Getting the results ---------------------------------------------- */
	double minDistance_ff = INFINITY;
	double maxDistance_ff = 0;
	double minDistance    = INFINITY;
	double maxDistance    = 0;
	int minDistanceIdx_ff, maxDistanceIdx_ff;
       	int minDistanceIdx,    maxDistanceIdx;

	for(int i = 0; i < ProfileLength; i++)
	{
		if(ff_get_double(&profile_ff[i]) < minDistance_ff && 
				ff_get_double(&profile_ff[i]) > 0)
		{
			minDistance_ff    = ff_get_double(&profile_ff[i]);
			minDistanceIdx_ff = profileIdxs_ff[i];
		}

		if(ff_get_double(&profile_ff[i]) > maxDistance_ff 
				&& ff_get_double(&profile_ff[i]) > 0)
		{
			maxDistance_ff    = ff_get_double(&profile_ff[i]);
			maxDistanceIdx_ff = profileIdxs_ff[i];
		}
		if(profile[i] < minDistance && profile[i] > 0)
		{
			minDistance    = profile[i];
			minDistanceIdx = profileIdxs[i];
		}

		if(profile[i] > maxDistance && profile[i] > 0)
		{
			maxDistance    = profile[i];
			maxDistanceIdx = profileIdxs[i];
		}
	}

	minDistance_ff = sqrt(minDistance_ff);
	maxDistance_ff = sqrt(maxDistance_ff);
	minDistance = sqrt(minDistance);
	maxDistance = sqrt(maxDistance);
	printf("[INFO] Results:\n");
	printf("----------------------------------------------\n");
	printf(" SCRIMP FF Min: %f Idx: %d\n", minDistance_ff, 
			minDistanceIdx_ff);
	printf(" SCRIMP FF Max: %f Idx: %d\n", maxDistance_ff, 
			maxDistanceIdx_ff);
	printf(" SCRIMP    Min: %f Idx: %d\n", minDistance, minDistanceIdx);
	printf(" SCRIMP    Max: %f Idx: %d\n", maxDistance, maxDistanceIdx);
	printf("----------------------------------------------\n");
	printf("[INFO] FlexFloat stats:\n");
	printf("----------------------------------------------\n");
	ff_print_stats();
	printf("----------------------------------------------\n");
	/* ------------------------------------------------------------------ */

	/* Saving results to file ------------------------------------------- */
	path_result = malloc(1000 * sizeof(char));
	strcpy(path_result, "result_");
	strcat(path_result, argv[1]);
	path_result[strlen(path_result)-4] = 0;
	strcat(path_result, ".csv");
	printf("[INFO] Saving %s ...\n", path_result);

	strcpy(path_result, PATH_RESULT);
	strcat(path_result, argv[1]);
	path_result[strlen(path_result)-4] = 0;
	strcat(path_result, ".csv"); 

	fp = fopen(path_result,"w");
	double error;
	for (int i = 0; i < ProfileLength; i++) {
		error = (fabs((sqrt(profile[i])
                                - sqrt(ff_get_double(&profile_ff[i]))))
                                / sqrt(profile[i])) * 100;
		fprintf(fp, "%d,%f,%f,%d,%f,%d,%f\n", i, tSeries[i], 
				sqrt(ff_get_double(&profile_ff[i])), 
				profileIdxs_ff[i], sqrt(profile[i]),
				profileIdxs[i], error);
	}
	fclose(fp);

	printf("[INFO] DONE\n");
	/* ------------------------------------------------------------------ */

	free(path_tSeries);
	free(path_result);
	free(tSeries);
	free(AMean); 
	free(ASigma);
	free(profile);
	free(profileIdxs);
	free(profileIdxs_ff);
	free(idx);
	free(tSeries_ff);
	free(AMean_ff);
	free(ASigma_ff);
	free(profile_ff);

	printf("##############################################\n");

	return 0;
}
