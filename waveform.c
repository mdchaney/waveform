/*
	waveform

	By Michael Chaney
	Copyright 2014 Michael Chaney Consulting Corporation - All rights reserved

	Released under terms of the MIT license or GNU General Public
	License v2.  See MIT-LICENSE or GPLv2-LICENSE files for details.

	The "points" option tells how many data points will be computed in
	total.  Bresenham's algorithm will be used to determine the number
	of samples used to compute each data point.  An RMS calculation will
	be used to compute each sample point, however the "--mean" flag can
	be used to cause a standard mean to be used.

	This program will take either a filename argument or read from stdin.
	As long as the wav/aif is in correct order (COMM/fmt chunk before
	SSND/data chunk) it'll work fine.

	I don't use libsndfile or other such libraries to make reading
	easier.  Instead, this is built to be as fast as possible and thus
	has separate code paths depending on what is needed to read the file.

	There are three optimizations here:

	1. Inner loops are as simple as possible with minimal branching
	inside of loops.

	2. Inlining of commonly used functions, and minimal use of function
	calls.

	3. Integer arithmetic where possible.  The RMS function does make use
	of some 64-bit integers but only uses floating point for the few
	operations where it is required.  This will not run as fast on a
	32-bit architecture and should probably be redone if you wish to use
	it on such a machine.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <math.h>

/* My structures are representation of file data structures and must be
 * used exactly as shown. */
#pragma pack(1)

static int verbose_flag = 0, debug_flag = 0, mono_flag = 0, raw_flag = 0;
static int help_flag = 0;
static unsigned long points = 1000;
static double scale = 256.0;
static int use_peak = 0, use_mean = 0, use_rms = 0;

typedef enum { SIGNED, UNSIGNED } Signing_t;
typedef enum { BIG, LITTLE, MIXED } Endianness_t;
typedef enum { RIFF_FILE, FORM_FILE } FileFormat_t;
typedef enum { WAVE_FILE, AIFF_FILE, AIFC_FILE } AudioFormat_t;
typedef enum { RMS, MEAN, PEAK } Algo_t;

/*
	These inline functions take care of byte swapping if the endiannesses
	of the file and machine differ.
*/

/* byte swap unsigned 16-bit int */
static inline uint16_t swap_uint16(uint16_t val) {
	return (val << 8) | (val >> 8 );
}

/* byte swap signed 16-bit int */
static inline int16_t swap_int16(int16_t val) {
	return (val << 8) | ((val >> 8) & 0xFF);
}

/* byte swap uint 32 */
static inline uint32_t swap_uint32(uint32_t val) {
	val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF ); 
	return (val << 16) | (val >> 16);
}

/* byte swap int 32 */
static inline int32_t swap_int32(int32_t val) {
	val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
	return (val << 16) | ((val >> 16) & 0xFFFF);
}

/* macro to handle 24-bit integer - PTR is the pointer to raw data, TARGET
   is an int32_t where we'll store the data.  PTR is assumed to be (uint8_t *).
	Also assuming that the bit shift operator will drag the sign bit along.
	The first two are for little-endian architectures, second two are for
	big-endian.
*/
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define get_unswapped_24_bit_int(PTR, TARGET) (TARGET=0, ((uint8_t *)&TARGET)[1] = *(PTR++), ((uint8_t *)&TARGET)[2] = *(PTR++), ((uint8_t *)&TARGET)[3] = *(PTR++), (TARGET = TARGET >> 8))

#define get_swapped_24_bit_int(PTR, TARGET) (TARGET=0, ((uint8_t *)&TARGET)[3] = *(PTR++), ((uint8_t *)&TARGET)[2] = *(PTR++), ((uint8_t *)&TARGET)[1] = *(PTR++), (TARGET = TARGET >> 8))

#else
#define get_unswapped_24_bit_int(PTR, TARGET) (TARGET=0, ((uint8_t *)&TARGET)[0] = *(PTR++), ((uint8_t *)&TARGET)[1] = *(PTR++), ((uint8_t *)&TARGET)[2] = *(PTR++), (TARGET = TARGET >> 8))

#define get_swapped_24_bit_int(PTR, TARGET) (TARGET=0, ((uint8_t *)&TARGET)[2] = *(PTR++), ((uint8_t *)&TARGET)[1] = *(PTR++), ((uint8_t *)&TARGET)[0] = *(PTR++), (TARGET = TARGET >> 8))
#endif

/* Here's where we implement Bresenham's algorithm for line
 * drawing.  Let's assume a line in the first octant (in polar
 * coordinates the angle is 0-45 degrees) where "x" is the
 * actual sample_count and "y" is the "points" variable.  We
 * assume then that "sample_count > points".
 *
 * In standard Bresenham talk that means that "sample_count" is
 * "dx" and "points" is "dy".  Using integer division the
 * number of samples per point will vary between two numbers,
 * say "lower_points_per_sample" and
 * "lower_points_per_sample+1".  When to use
 * "lower_points_per_sample" and when to use
 * "lower_points_per_sample+1" is the trick.
 *
 * Let's say there are 9000 samples and we want 1000 points.
 * In that case each point will be computed from 9 samples.
 *
 * If we add 1 sample then we're at 9001 samples.  In that case
 * we will need 999 groups of 9 samples and one group of 10
 * samples.  If we have 9005 samples then it'll be 995 groups
 * of 9 samples and 5 groups of 10 samples.
 *
 *    Samples            9 samples/group    10 samples/group
 *       9000                  1000                   0
 *       9001                   999                   1
 *       9005                   995                   5
 *       9499                   501                 499
 *       9500                   500                 500
 *       9999                     1                 999
 *
 * Here's where it gets weird.  We can actually again employ
 * Bresenham's algorithm to determine the interval at which we
 * use "lower_points_per_sample" samples and when we use
 * "lower_points_per_sample+1".  As the "slope" exceeds 45
 * degrees (at "points / 2") we have to use
 * "lower_points_per_sample+1" more often than
 * "lower_points_per_sample".
 *
 * Take the first case, where samples is between 9000 and 9499.
 * Generalized (with integer arithmetic):
 *
 * lower_points_per_sample = samples / points leftover =
 * samples - (lower_points_per_sample * points)
 *
 * At this point with the example numbers above "samples" is
 * "1000" * and "leftover" is the last column.  For our example
 * we examine numbers where "leftover / points < 0.5" for
 * simplicity.  In that case we need to stick between 0 and 499
 * "lower_points_per_sample+1" sample groups in with the other
 * samples.  So every so often we have to grab
 * "lower_points_per_sample+1" instead of
 * "lower_points_per_sample".
 *
 */

int* get_sample_group_sizes(int sample_count, int points) {

	long lower_points_per_sample = sample_count / points;
	long leftover = sample_count - (lower_points_per_sample * points) + 1;

	if (debug_flag) fprintf(stderr, "lower_points_per_sample: %ld, leftover: %ld\n", lower_points_per_sample, leftover);

	long samples_left = sample_count;
	int *sample_group_sizes;

	sample_group_sizes = (int *) malloc(sizeof(int) * (points+1));

	int i = 0, j = 0, k = 0, jump_counter;

	jump_counter = leftover - points;

	while (samples_left > 0 && i < points) {
		if (jump_counter > 0) {
			samples_left -= lower_points_per_sample+1;
			sample_group_sizes[i] = lower_points_per_sample+1;
			jump_counter -= points;
		} else {
			samples_left -= lower_points_per_sample;
			sample_group_sizes[i] = lower_points_per_sample;
		}
		i++;
		jump_counter += leftover;
	}

	if (debug_flag) {
		samples_left = sample_count;
		for (i=0 ; i<points; i++) {
			fprintf(stderr, "%5d:   %10d   %10ld\n", i, sample_group_sizes[i], samples_left);
			samples_left -= sample_group_sizes[i];
		}
		fprintf(stderr, "%5d:                %10ld\n", i, samples_left);
	}

	return(sample_group_sizes);
}

/*
	This group of functions handles 8-bit samples.  They can be:
	1. peak, mean, or rms
	2. 1 or 2 channel
	3. signed or unsigned (WAV uses unsigned ints offset by 128)
	There are 6 functions in all, as signed/unsigned functions have been
	rolled in together.
*/

int waveform_1_channel_8_bit_peak(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 1
	 * bits_per_sample == 8
	 * algorithm == PEAK
	 */

	int j;

	int16_t sample_point, peak_0=0;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t*)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs((int16_t)*u_sample_pointer - 128);
			if (sample_point > peak_0) peak_0 = sample_point;
			u_sample_pointer++;
		}
	}

	printf("%u\n", (int)floor((double)peak_0 * scale / 128.0));

	return(1);
}

int waveform_2_channel_8_bit_peak(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 2
	 * bits_per_sample == 8
	 * algorithm == PEAK
	 */

	int j;

	int16_t sample_point, peak_0=0, peak_1=0;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_1) peak_1 = sample_point;
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t *)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs((int16_t)*u_sample_pointer - 128);
			if (sample_point > peak_0) peak_0 = sample_point;
			u_sample_pointer++;
			sample_point = abs((int16_t)*u_sample_pointer - 128);
			if (sample_point > peak_1) peak_1 = sample_point;
			u_sample_pointer++;
		}
	}

	if (mono_flag) {
		if (peak_1 > peak_0) { peak_0 = peak_1; }
		printf("%u\n", (int)floor((double)peak_0 * scale / 128.0));
	} else {
		printf("%u,%u\n", (int)floor((double)peak_0 * scale / 128.0), (int)floor((double)peak_1 * scale / 128.0));
	}

	return(1);
}

int waveform_1_channel_8_bit_mean(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 1
	 * bits_per_sample == 8
	 * algorithm == MEAN
	 */

	int j;

	int64_t sum_of_samples_0=0;
	double mean_0;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs((int64_t)*sample_pointer);
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t*)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs((int64_t)*u_sample_pointer - 128);
			u_sample_pointer++;
		}
	}

	mean_0 = (double)sum_of_samples_0 / (double)sample_group_size;
	printf("%u\n", (int)floor((double)mean_0 * scale / 128.0));

	return(1);
}

int waveform_2_channel_8_bit_mean(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 2
	 * bits_per_sample == 8
	 * algorithm == MEAN
	 */

	int j;

	int64_t sample_point, sum_of_samples_0=0, sum_of_samples_1=0;
	double mean_0, mean_1;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs((int64_t)*sample_pointer);
			sample_pointer++;
			sum_of_samples_1 += abs((int64_t)*sample_pointer);
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t*)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs((int64_t)*u_sample_pointer - 128);
			u_sample_pointer++;
			sum_of_samples_1 += abs((int64_t)*u_sample_pointer - 128);
			u_sample_pointer++;
		}
	}

	if (mono_flag) {
		mean_0 = (double)(sum_of_samples_0 + sum_of_samples_1) / ((double)sample_group_size * 2.0);
		printf("%u\n", (int)floor((double)mean_0 * scale / 128.0));
	} else {
		mean_0 = (double)sum_of_samples_0 / (double)sample_group_size;
		mean_1 = (double)sum_of_samples_1 / (double)sample_group_size;
		printf("%u,%u\n", (int)floor((double)mean_0 * scale / 128.0), (int)floor((double)mean_1 * scale / 128.0));
	}

	return(1);
}

int waveform_1_channel_8_bit_rms(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 1
	 * bits_per_sample == 8
	 * algorithm == RMS
	 */

	int j;

	int64_t sample_point_64, sum_of_squares_0=0;
	double rms_0;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point_64 = (int64_t)*sample_pointer;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t*)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point_64 = (int64_t)*u_sample_pointer - 128;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			u_sample_pointer++;
		}
	}

	rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);
	printf("%u\n", (int)floor((double)rms_0 * scale / 128.0));

	return(1);
}

int waveform_2_channel_8_bit_rms(int8_t *samples, int sample_group_size, Signing_t file_signing) {

	/* this code is where:
	 * data format == signed/unsigned (see file_signing)
	 * channel_count == 2
	 * bits_per_sample == 8
	 * algorithm == RMS
	 */

	int j;

	int64_t sample_point_64, sum_of_squares_0=0, sum_of_squares_1=0;
	double rms_0, rms_1;

	if (file_signing == SIGNED) {
		int8_t *sample_pointer = samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point_64 = (int64_t)*sample_pointer;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			sample_pointer++;
			sample_point_64 = (int64_t)*sample_pointer;
			sum_of_squares_1 += (sample_point_64 * sample_point_64);
			sample_pointer++;
		}
	} else {
		uint8_t *u_sample_pointer = (uint8_t *)samples;
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point_64 = (int64_t)*u_sample_pointer - 128;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			u_sample_pointer++;
			sample_point_64 = (int64_t)*u_sample_pointer - 128;
			sum_of_squares_1 += (sample_point_64 * sample_point_64);
			u_sample_pointer++;
		}
	}

	if (mono_flag) {
		rms_0 = sqrt((double)(sum_of_squares_0 + sum_of_squares_1) / ((double)sample_group_size * 2.0));
		printf("%u\n", (int)floor((double)rms_0 * scale / 128.0));
	} else {
		rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);
		rms_1 = sqrt((double)sum_of_squares_1 / (double)sample_group_size);
		printf("%u,%u\n", (int)floor((double)rms_0 * scale / 128.0), (int)floor((double)rms_1 * scale / 128.0));
	}

	return(1);
}

/*
	This group of functions handles 16-bit samples.  They can be:
	1. peak, mean, or rms
	2. 1 or 2 channel
	So, for 16-bit samples there are 6 functions, with different branches
	to handle endianness issues.
*/

int waveform_1_channel_16_bit_peak(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample == 16
	 * algorithm == PEAK
	 */

	int j;

	int16_t *sample_pointer = samples;

	int16_t sample_point, peak_0=0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(swap_int16(*sample_pointer));
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
		}
	}

	printf("%u\n", (int)floor((double)peak_0 * scale / 32768.0));

	return(1);
}

int waveform_2_channel_16_bit_peak(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample == 16
	 * algorithm == PEAK
	 */

	int j;

	int16_t *sample_pointer = samples;

	int16_t sample_point, peak_0=0, peak_1=0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
			sample_point = abs(*sample_pointer);
			if (sample_point > peak_1) peak_1 = sample_point;
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = abs(swap_int16(*sample_pointer));
			if (sample_point > peak_0) peak_0 = sample_point;
			sample_pointer++;
			sample_point = abs(swap_int16(*sample_pointer));
			if (sample_point > peak_1) peak_1 = sample_point;
			sample_pointer++;
		}
	}

	if (mono_flag) {
		if (peak_1 > peak_0) { peak_0 = peak_1; }
		printf("%u\n", (int)floor((double)peak_0 * scale / 32768.0));
	} else {
		printf("%u,%u\n", (int)floor((double)peak_0 * scale / 32768.0), (int)floor((double)peak_1 * scale / 32768.0));
	}

	return(1);
}

int waveform_1_channel_16_bit_mean(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample == 16
	 * algorithm == MEAN
	 */

	int j;

	int16_t *sample_pointer = samples;

	int64_t sample_point, sum_of_samples_0=0;
	double mean_0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs(*sample_pointer);
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs(swap_int16(*sample_pointer));
			sample_pointer++;
		}
	}

	mean_0 = (double)sum_of_samples_0 / ((double)sample_group_size * 2.0);

	printf("%u\n", (int)floor(mean_0 * scale / 32768.0));

	return(1);
}

int waveform_2_channel_16_bit_mean(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample == 16
	 * algorithm == MEAN
	 */

	int j;

	int16_t *sample_pointer = samples;

	int64_t sample_point, sum_of_samples_0=0, sum_of_samples_1=0;
	double mean_0, mean_1;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs(*sample_pointer);
			sample_pointer++;
			sum_of_samples_1 += abs(*sample_pointer);
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			sum_of_samples_0 += abs(swap_int16(*sample_pointer));
			sample_pointer++;
			sum_of_samples_1 += abs(swap_int16(*sample_pointer));
			sample_pointer++;
		}
	}

	if (mono_flag) {
		mean_0 = (double)(sum_of_samples_0 + sum_of_samples_1) / ((double)sample_group_size * 2.0);
		printf("%u\n", (int)floor(mean_0 * scale / 32768.0));
	} else {
		mean_0 = (double)sum_of_samples_0 / (double)sample_group_size;
		mean_1 = (double)sum_of_samples_1 / (double)sample_group_size;
		printf("%u,%u\n", (int)floor(mean_0 * scale / 32768.0), (int)floor(mean_1 * scale / 32768.0));
	}

	return(1);
}

int waveform_1_channel_16_bit_rms(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample == 16
	 * algorithm == RMS
	 */

	int j;
	int16_t *sample_pointer=samples;

	int64_t sample_point, sum_of_squares_0=0;
	double rms_0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = *sample_pointer;
			sum_of_squares_0 += (sample_point * sample_point);
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			sample_point = swap_int16(*sample_pointer);
			sum_of_squares_0 += (sample_point * sample_point);
			sample_pointer++;
		}
	}

	/* At this point we have the sums of the squares of
	 * the sample group.  We'll do our floating point math
	 * here to get the square root. */
	rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);

	printf("%u\n", (int)floor(rms_0 * scale / 32768.0));

	return(1);
}

int waveform_2_channel_16_bit_rms(int16_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample == 16
	 * algorithm == RMS
	 */

	int j;

	int16_t *sample_pointer = samples;

	int64_t sum_of_squares_0=0, sum_of_squares_1=0;
	double rms_0, rms_1;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			int64_t sample_point = *sample_pointer;
			sum_of_squares_0 += (sample_point * sample_point);
			sample_pointer++;
			sample_point = *sample_pointer;
			sum_of_squares_1 += (sample_point * sample_point);
			sample_pointer++;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			int64_t sample_point = swap_int16(*sample_pointer);
			sum_of_squares_0 += (sample_point * sample_point);
			sample_pointer++;
			sample_point = swap_int16(*sample_pointer);
			sum_of_squares_1 += (sample_point * sample_point);
			sample_pointer++;
		}
	}

	/* At this point we have the sums of the squares of
	 * the sample group.  We'll do our floating point math
	 * here to get the square root. */

	if (mono_flag) {
		rms_0 = sqrt((double)(sum_of_squares_0 + sum_of_squares_1) / ((double)sample_group_size * 2.0));
		printf("%u\n", (int)floor(rms_0 * scale / 32768.0));
	} else {
		rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);
		rms_1 = sqrt((double)sum_of_squares_1 / (double)sample_group_size);
		printf("%u,%u\n", (int)floor(rms_0 * scale / 32768.0), (int)floor(rms_1 * scale / 32768.0));
	}

	return(1);
}

/*
	These are the functions for 24-bit samples.
	1. peak, mean, or rms
	2. 1 or 2 channel
	So, for 24-bit samples there are 6 functions, with different branches
	to handle endianness issues.
*/

int waveform_1_channel_24_bit_peak(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample ==24 
	 * algorithm == PEAK
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point_0, peak_0=0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point_0);
			sample_point_0 = abs(sample_point_0);
			if (sample_point_0 > peak_0) peak_0 = sample_point_0;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point_0);
			sample_point_0 = abs(sample_point_0);
			if (sample_point_0 > peak_0) peak_0 = sample_point_0;
		}
	}
	printf("%u\n", (int)floor((double)peak_0 * scale / 8388608.0));

	return(1);
}

int waveform_2_channel_24_bit_peak(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample ==24 
	 * algorithm == PEAK
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point_0, sample_point_1, peak_0=0, peak_1=0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point_0);
			get_unswapped_24_bit_int(sample_pointer, sample_point_1);

			sample_point_0 = abs(sample_point_0);
			sample_point_1 = abs(sample_point_1);

			if (sample_point_0 > peak_0) peak_0 = sample_point_0;
			if (sample_point_1 > peak_1) peak_1 = sample_point_1;
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point_0);
			get_swapped_24_bit_int(sample_pointer, sample_point_1);

			sample_point_0 = abs(sample_point_0);
			sample_point_1 = abs(sample_point_1);

			if (sample_point_0 > peak_0) peak_0 = sample_point_0;
			if (sample_point_1 > peak_1) peak_1 = sample_point_1;
		}
	}

	if (mono_flag) {
		if (peak_1 > peak_0) peak_0 = peak_1;
		printf("%u\n", (int)floor((double)peak_0 * scale / 8388608.0));
	} else {
		printf("%u,%u\n", (int)floor((double)peak_0 * scale / 8388608.0), (int)floor((double)peak_1 * scale / 8388608.0));
	}

	return(1);
}

int waveform_1_channel_24_bit_mean(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample ==24 
	 * algorithm == MEAN
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point_0;
	int64_t sum_of_samples_0=0;
	double mean_0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point_0);
			sum_of_samples_0 += abs(sample_point_0);
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point_0);
			sum_of_samples_0 += abs(sample_point_0);
		}
	}

	mean_0 = (double)sum_of_samples_0 / (double)sample_group_size;
	printf("%u\n", (int)floor(mean_0 * scale / 8388608.0));

	return(1);
}

int waveform_2_channel_24_bit_mean(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample ==24 
	 * algorithm == MEAN
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point_0, sample_point_1;
	int64_t sum_of_samples_0=0, sum_of_samples_1=0;
	double mean_0, mean_1;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point_0);
			get_unswapped_24_bit_int(sample_pointer, sample_point_1);
			sum_of_samples_0 += abs(sample_point_0);
			sum_of_samples_1 += abs(sample_point_1);
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point_0);
			get_swapped_24_bit_int(sample_pointer, sample_point_1);
			sum_of_samples_0 += abs(sample_point_0);
			sum_of_samples_1 += abs(sample_point_1);
		}
	}

	if (mono_flag) {
		mean_0 = (double)(sum_of_samples_0 + sum_of_samples_1) / ((double)sample_group_size * 2.0);
		printf("%u\n", (int)floor(mean_0 * scale / 8388608.0));
	} else {
		mean_0 = (double)sum_of_samples_0 / (double)sample_group_size;
		mean_1 = (double)sum_of_samples_1 / (double)sample_group_size;
		printf("%u,%u\n", (int)floor(mean_0 * scale / 8388608.0), (int)floor(mean_1 * scale / 8388608.0));
	}

	return(1);
}

int waveform_1_channel_24_bit_rms(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 1
	 * bits_per_sample ==24 
	 * algorithm == RMS
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point;
	int64_t sample_point_64;
	int64_t sum_of_squares_0=0;
	double rms_0;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
		}
	}

	/* At this point we have the sums of the squares of
	 * the sample group.  We'll do our floating point math
	 * here to get the square root. */

	rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);
	printf("%u\n", (int)floor(rms_0 * scale / 8388608.0));

	return(1);
}

int waveform_2_channel_24_bit_rms(uint8_t *samples, int sample_group_size, int same_endianness) {

	/* this code is where:
	 * channel_count == 2
	 * bits_per_sample ==24 
	 * algorithm == RMS
	 */

	/* Note that 8388608 is 2^23 */

	int j;

	uint8_t *sample_pointer = samples;

	int32_t sample_point;
	int64_t sample_point_64;
	int64_t sum_of_squares_0=0, sum_of_squares_1=0;
	double rms_0, rms_1;

	if (same_endianness) {
		for (j=0 ; j<sample_group_size ; j++) {
			get_unswapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			get_unswapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_1 += (sample_point_64 * sample_point_64);
		}
	} else {
		for (j=0 ; j<sample_group_size ; j++) {
			get_swapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_0 += (sample_point_64 * sample_point_64);
			get_swapped_24_bit_int(sample_pointer, sample_point);
			sample_point_64 = sample_point;
			sum_of_squares_1 += (sample_point_64 * sample_point_64);
		}
	}

	/* At this point we have the sums of the squares of
	 * the sample group.  We'll do our floating point math
	 * here to get the square root. */

	if (mono_flag) {
		rms_0 = sqrt((double)(sum_of_squares_0 + sum_of_squares_1) / ((double)sample_group_size * 2.0));
		printf("%u\n", (int)floor(rms_0 * scale / 8388608.0));
	} else {
		rms_0 = sqrt((double)sum_of_squares_0 / (double)sample_group_size);
		rms_1 = sqrt((double)sum_of_squares_1 / (double)sample_group_size);
		printf("%u,%u\n", (int)floor(rms_0 * scale / 8388608.0), (int)floor(rms_1 * scale / 8388608.0));
	}

	return(1);
}

/*
	The following functions are dispatchers and are based on number of
	channels and sample size.  This leaves us with 6 functions since we
	do 8, 16, and 24 bit samples and one or two channels.
*/

int waveform_1_channel_8_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Signing_t file_signing) {
	int i, items_read;

	int (*funcptr)(int8_t*, int, Signing_t);

	if (algorithm == RMS) {
		funcptr = &waveform_1_channel_8_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_1_channel_8_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_1_channel_8_bit_mean;
	}

	int8_t *samples;
	samples = (int8_t *) malloc(sample_group_sizes[0]+1);  /* 1 channel, 1 byte samples */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 1, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, file_signing);
	}

	free(samples);

	return(1);
}

int waveform_1_channel_16_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness) {

	int i, items_read;

	int (*funcptr)(int16_t*, int, int);

	if (algorithm == RMS) {
		funcptr = &waveform_1_channel_16_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_1_channel_16_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_1_channel_16_bit_mean;
	}

	int16_t *samples;
	samples = (int16_t *) malloc(2 * (sample_group_sizes[0]+1));  /* 1 channel, 2 byte samples */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 2, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, (machine_endianness == data_endianness ? 1 : 0));
	}

	free(samples);

	return(1);
}

int waveform_1_channel_24_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness) {

	int i, items_read;

	int (*funcptr)(uint8_t*, int, int);

	if (algorithm == RMS) {
		funcptr = &waveform_1_channel_24_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_1_channel_24_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_1_channel_24_bit_mean;
	}

	uint8_t *samples;
	samples = (uint8_t *) malloc(3 * (sample_group_sizes[0]+1));  /* 1 channel, 3 byte samples */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 3, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, (machine_endianness == data_endianness ? 1 : 0));
	}

	free(samples);

	return(1);
}

int waveform_2_channel_8_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Signing_t file_signing) {

	int i, items_read;

	int (*funcptr)(int8_t*, int, Signing_t);

	if (algorithm == RMS) {
		funcptr = &waveform_2_channel_8_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_2_channel_8_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_2_channel_8_bit_mean;
	}

	int8_t *samples;
	samples = (int8_t *) malloc(2 * (sample_group_sizes[0]+1));  /* 2 channels, 1 byte samples */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 2, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, file_signing);
	}

	free(samples);

	return(1);
}

int waveform_2_channel_16_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness) {

	int i, items_read;

	int (*funcptr)(int16_t*, int, int);

	if (algorithm == RMS) {
		funcptr = &waveform_2_channel_16_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_2_channel_16_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_2_channel_16_bit_mean;
	}

	int16_t *samples;
	samples = (int16_t *) malloc(2 * 2 * (sample_group_sizes[0]+1));  /* 2 channels, 2 byte samples */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 4, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, (machine_endianness == data_endianness ? 1 : 0));
	}

	free(samples);

	return(1);
}

int waveform_2_channel_24_bit(FILE *fd, int *sample_group_sizes, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness) {

	int i, items_read;

	int (*funcptr)(uint8_t*, int, int);

	if (algorithm == RMS) {
		funcptr = &waveform_2_channel_24_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_2_channel_24_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_2_channel_24_bit_mean;
	}

	uint8_t *samples;
	samples = (uint8_t *) malloc(2 * 3 * (sample_group_sizes[0]+16));  /* 2 channels, 3 byte samples, plus padding */

	for (i=0 ; i<points; i++) {

		int points_to_read = sample_group_sizes[i];

		items_read = fread(samples, 2 * 3, points_to_read, fd);

		if (items_read != points_to_read) {
			if (feof(fd)) {
				fprintf(stderr, "Unexpected EOF\n");
				exit(1);
			} else {
				fprintf(stderr, "Chunk unreadable\n");
				exit(1);
			}
		}

		(*funcptr)(samples, points_to_read, (machine_endianness == data_endianness ? 1 : 0));
	}

	free(samples);

	return(1);
}

/*
  These are the second level dispatchers which are based only on
  sample size.
*/

int waveform_1_channel(FILE *fd, int *sample_group_sizes, int bits_per_sample, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness, Signing_t file_signing) {
	if (bits_per_sample == 8) {
		waveform_1_channel_8_bit(fd, sample_group_sizes, algorithm, file_signing);
	} else if (bits_per_sample == 16) {
		waveform_1_channel_16_bit(fd, sample_group_sizes, algorithm, machine_endianness, data_endianness);
	} else if (bits_per_sample == 24) {
		waveform_1_channel_24_bit(fd, sample_group_sizes, algorithm, machine_endianness, data_endianness);
	}
	return 0;
}

int waveform_2_channel(FILE *fd, int *sample_group_sizes, int bits_per_sample, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness, Signing_t file_signing) {
	if (bits_per_sample == 8) {
		waveform_2_channel_8_bit(fd, sample_group_sizes, algorithm, file_signing);
	} else if (bits_per_sample == 16) {
		waveform_2_channel_16_bit(fd, sample_group_sizes, algorithm, machine_endianness, data_endianness);
	} else if (bits_per_sample == 24) {
		waveform_2_channel_24_bit(fd, sample_group_sizes, algorithm, machine_endianness, data_endianness);
	}
	return 0;
}

/*
  This is the top level dispatcher.  It will dispatch based only on
  channel count.
*/
int calculate_waveform(FILE *fd, int sample_count, int channel_count, int bits_per_sample, Algo_t algorithm, Endianness_t machine_endianness, Endianness_t data_endianness, Signing_t file_signing) {

	int* sample_group_sizes = get_sample_group_sizes(sample_count, points);

	if (channel_count == 1) {
		return waveform_1_channel(fd, sample_group_sizes, bits_per_sample, algorithm, machine_endianness, data_endianness, file_signing);
	} else if (channel_count == 2) {
		return waveform_2_channel(fd, sample_group_sizes, bits_per_sample, algorithm, machine_endianness, data_endianness, file_signing);
	}

	free(sample_group_sizes);

	return(0);
}

/* This is confusing because an item in a chunk is a sample pair, and
 * each sample is 2 bytes. */

int calculate_waveform_from_raw(int16_t **chunks, int items_in_chunk, int sample_count, Algo_t algorithm) {
	int* sample_group_sizes = get_sample_group_sizes(sample_count, points);

	int i, items_read;

	int (*funcptr)(int16_t*, int, int);

	if (algorithm == RMS) {
		funcptr = &waveform_2_channel_16_bit_rms;
	} else if (algorithm == PEAK) {
		funcptr = &waveform_2_channel_16_bit_peak;
	} else if (algorithm == MEAN) {
		funcptr = &waveform_2_channel_16_bit_mean;
	}

	int16_t *samples;
	samples = (int16_t *) malloc(2 * 2 * (sample_group_sizes[0]+1));  /* 2 channels, 2 byte samples */

	int points_to_read, samples_offset = 0, current_chunk = 0, chunk_offset = 0;

	for (i=0 ; i<points; i++) {

		points_to_read = sample_group_sizes[i];

		/* I know that points_to_read will always be less than the size of
		 * a chunk. */

		if (debug_flag) fprintf(stderr, "%5d  %8d %9d  %3d  %7d\n", i, points_to_read, samples_offset, current_chunk, chunk_offset);

		if (items_in_chunk - chunk_offset < points_to_read) {
			/* need to pull from this chunk and next */
			memcpy(samples, &chunks[current_chunk][chunk_offset*2], (items_in_chunk - chunk_offset) * 4);
			if (debug_flag) fprintf(stderr, "split  %7d\n", items_in_chunk - chunk_offset);
			current_chunk++;
			memcpy(&samples[(items_in_chunk - chunk_offset)*2], &chunks[current_chunk][0], (points_to_read - (items_in_chunk - chunk_offset)) * 4);
			if (debug_flag) fprintf(stderr, "split  %7d\n", points_to_read - (items_in_chunk - chunk_offset));
			chunk_offset = points_to_read - (items_in_chunk - chunk_offset);
			(*funcptr)(samples, points_to_read, 1);
		} else {
			(*funcptr)(&chunks[current_chunk][chunk_offset*2], points_to_read, 1);
			chunk_offset += points_to_read;
			if (chunk_offset >= items_in_chunk) {
				chunk_offset = 0;
				current_chunk++;
			}
		}

		samples_offset += points_to_read;
	}

	free(samples);

	free(sample_group_sizes);

	return(1);
}

int main(int argc, char **argv) {
	int i = 0, j = 0, k = 0; /* for use in loops */

	int c;

	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 },
			{ "debug", no_argument, &debug_flag, 1 },
			{ "mono", no_argument, &mono_flag, 1 },
			{ "peak", no_argument, &use_peak, 1 },
			{ "mean", no_argument, &use_mean, 1 },
			{ "rms", no_argument, &use_rms, 1 },
			{ "points", required_argument, 0, 'p' },
			{ "scale", required_argument, 0, 's' },
			{ "raw", no_argument, &raw_flag, 1 },
			{ "help", no_argument, &help_flag, 1 },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		c = getopt_long(argc, argv, "", long_options, &option_index);

		if (c == -1) break;

		switch(c) {
			case 0:
				if (long_options[option_index].flag != 0) break;
				break;

			case 'p':
				/* parse optarg to get # */
				if (sscanf(optarg, "%lu", &points) != 1) {
					fprintf(stderr, "Expecting number for 'points', got `%s'\n", optarg);
					exit(1);
				}
				break;

			case 's':
				/* parse optarg to get # */
				if (sscanf(optarg, "%lf", &scale) != 1) {
					fprintf(stderr, "Expecting number for 'scale', got `%s'\n", optarg);
					exit(1);
				}
				break;

			case '?': /* error condition */
				break;

			default:
				abort();
		}
	}

	if (help_flag) {
		fprintf(stdout,
"Usage: %s [--verbose] [--mono] [--peak|--mean|--rms] [--points x] [--scale y] [--raw] [filename]\n\
\n\
--verbose          - show more verbose output\n\
--mono             - output only one channel, even if input is stereo\n\
--peak             - choose the peak value from each sample group\n\
--mean             - calculate and output the mean value from each sample group\n\
--rms              - calculate and output the RMS (root mean squared) value \n\
                     from each sample group\n\
--points x         - output x points (defaults to 1000)\n\
--scale y          - scale the output values from 0 to y, inclusive\n\
--raw              - read raw PCM data (output of `lame --decode -t ...')\n\
\n\
See README for more details.\n", argv[0]);
		exit(0);
	}

	if (!use_peak && !use_mean && !use_rms) use_rms = 1;

	Algo_t algorithm;

	if (use_peak) {
		algorithm = PEAK;
	} else if (use_mean) {
		algorithm = MEAN;
	} else if (use_rms) {
		algorithm = RMS;
	}

	if (verbose_flag) {
		if (debug_flag) fprintf(stderr, "verbose!\n");
	}

	if (points) {
		if (debug_flag) fprintf(stderr, "points: %lu\n", points);
	}

	if (algorithm == MEAN) {
		if (debug_flag) fprintf(stderr, "using mean\n");
	} else if (algorithm == RMS) {
		if (debug_flag) fprintf(stderr, "using rms\n");
	} else if (algorithm == PEAK) {
		if (debug_flag) fprintf(stderr, "using peak\n");
	} else {
		if (debug_flag) fprintf(stderr, "confused about algorithm\n");
	}

	char* filename;

	if (debug_flag) fprintf(stderr, "optind %d, argc %d\n", optind, argc);

	if (optind-1 < argc) {
		filename = argv[optind];
	}

	if (filename) {
		if (debug_flag) fprintf(stderr, "Filename: %s\n", filename);
	}

	if (!filename) {
		if (debug_flag) fprintf(stderr, "Reading stdin\n");
	}

	if (filename) {
		if (!freopen(filename, "r", stdin)) {
			fprintf(stderr, "Cannot open `%s' for reading\n", filename);
			exit(1);
		}
	}

	Signing_t file_signing;
	Endianness_t machine_endianness, file_endianness, data_endianness;

	uint32_t checker=0x01234567;

	switch (*((uint8_t*)(&checker))) {
		case 0x01:
			machine_endianness = BIG;
			break;
		case 0x23:
			machine_endianness = MIXED;
			break;
		default:
			machine_endianness = LITTLE;
			break;
	}

	if (debug_flag) fprintf(stderr, "Your architecture is ");

	if (machine_endianness == BIG) {
		if (debug_flag) fprintf(stderr, "big endian\n");
	} else if (machine_endianness == LITTLE) {
		if (debug_flag) fprintf(stderr, "little endian\n");
	} else if (machine_endianness == MIXED) {
		if (debug_flag) fprintf(stderr, "mixed endian\n");
	}

	size_t items_read;

	if (raw_flag) {

		/* reading raw pcm - probably from lame.  For now we'll assume
		 * 16-bit stereo with same endianness as the machine.  In this
		 * case we don't know how long the file will be, so we'll read it
		 * all in to memory so that we'll know how big it is, then we'll
		 * run the proper algorithm.  We'll read it in 1MB chunks.  */

		int64_t total_chunks = 0, item_size = 4, items_in_chunk = 250000, sample_count;

		/* Hard coding a maximum size here of 200 chunks.  That would be
		 * almost 20 minutes of audio, something of that size should be
		 * rare. */
		int chunk_limit=200;
		int16_t *chunks[chunk_limit];

		while (total_chunks<chunk_limit) {
			chunks[total_chunks] = malloc(item_size * items_in_chunk);

			items_read = fread(chunks[total_chunks], item_size, items_in_chunk, stdin);

			sample_count += items_read;
			total_chunks++;
	
	   	if (items_read < items_in_chunk) {
				if (feof(stdin)) {
					break;
				} else {
		      	fprintf(stderr, "Some sort of weird read error\n");
			      exit(1);
			   }
			}
		}

		if (total_chunks >= chunk_limit) {
			fprintf(stderr, "Too much data\n");
			exit(1);
		}

		/* We now have a bunch of samples, time to do whatever with them. */

		if (debug_flag) fprintf(stderr, "Sample Count: %d\n", sample_count);
	
		calculate_waveform_from_raw(chunks, items_in_chunk, sample_count, algorithm);

		for (i=0 ; i<total_chunks; i++) {
			free(chunks[i]);
		}
	
	} else {

		/* At this point, options are parsed and the file is opened.  We will
		 * read the first 8 bytes and see if it's a RIFF header. */
	
		typedef struct RIFF_HEADER {
			char RIFF[4];
			int32_t file_size_sans_riff_header;
		} riff_header_struct;
	
		riff_header_struct riff_header;
	
		items_read = fread(&riff_header, sizeof(riff_header_struct), 1, stdin);
	
		if (items_read != 1) {
			fprintf(stderr, "Header unreadable\n");
			exit(1);
		}
	
		if (strncmp(riff_header.RIFF, "RIFF", 4) != 0 && strncmp(riff_header.RIFF, "FORM", 4) != 0) {
			fprintf(stderr, "No RIFF/FORM header\n");
			exit(1);
		}
	
		FileFormat_t file_format;
	
		/* If it's a FORM then the file size will be in big-endian order */
		if (strncmp(riff_header.RIFF, "FORM", 4) == 0) {
			file_format = FORM_FILE;
		} else {
			file_format = RIFF_FILE;
		}
	
		if (file_format == FORM_FILE) {
			if (debug_flag) fprintf(stderr, "FORM file format\n");
			file_endianness = BIG;
			file_signing = SIGNED;
		} else if (file_format == RIFF_FILE) {
			if (debug_flag) fprintf(stderr, "RIFF file format\n");
			file_endianness = LITTLE;
			/* file signing depends on sample size */
		} else {
			if (debug_flag) fprintf(stderr, "Confusion in file format\n");
			exit(1);
		}
	
		if (debug_flag) fprintf(stderr, "Your file is ");
	
		if (file_endianness == BIG) {
			if (debug_flag) fprintf(stderr, "big endian\n");
		} else if (file_endianness == LITTLE) {
			if (debug_flag) fprintf(stderr, "little endian\n");
		}
	
		if (file_endianness != machine_endianness) {
			if (debug_flag) fprintf(stderr, "File size (sans RIFF/FORM header): %d\n", swap_int32(riff_header.file_size_sans_riff_header));
		} else {
			if (debug_flag) fprintf(stderr, "File size (sans RIFF/FORM header): %d\n", riff_header.file_size_sans_riff_header);
		}
	
		/* Now, we need to see what kind of file it is.  It's a RIFF if it
		 * got this far, this should be WAVE, AIFF or AIFC. */
	
		char audio_type[4];
	
		items_read = fread(audio_type, 4, 1, stdin);
	
		AudioFormat_t audio_format;
	
		if (debug_flag) fprintf(stderr, "File type: %c%c%c%c\n", audio_type[0], audio_type[1], audio_type[2], audio_type[3]);
	
		if (strncmp(audio_type, "WAVE", 4) == 0) {
			audio_format = WAVE_FILE;
			data_endianness = LITTLE;
		} else if (strncmp(audio_type, "AIFF", 4) == 0) {
			audio_format = AIFF_FILE;
			data_endianness = BIG;
		} else if (strncmp(audio_type, "AIFC", 4) == 0) {
			audio_format = AIFC_FILE;
			data_endianness = LITTLE;
		} else {
			fprintf(stderr, "Unknown file type %c%c%c%c", audio_type[0], audio_type[1], audio_type[2], audio_type[3]);
			exit(1);
		}
	
		/* Now, we have the file type, machine type, and audio type - time
			to read chunks.
	
		   We're only going to accept files for now where the description
		   precedes the data.  This is normal, anyway.
		 
		 	For AIFF/AIFC:
	
		 	COMM chunk - describes file
		 	SSND chunk - holds actual sample data
	
			For WAVE:
	
			fmt\0 chunk - describes file
			data chunk  - holds actual sample data
	
			For WAVE and AIFC (with "sowt" compression) the data will be in
			little-endian order.  For AIFF with no compression the data will
			be in big-endian order.
	
			If sample size is 8 bits in a WAVE file the samples are unsigned
			8-bit integers with 128 being "0".  To get signed amplitude we
			need to subtract 128 and turn it into a signed 8-bit integer
			(range of -128 through +127).
	
			If sample size is 16 bits the samples are already signed 16-bit
			integers and can be handled directly.
	
			For speed we have to simply handle each case separately.
		*/
	
		int we_are_done = 0;
	
		int bits_per_sample, channel_count, sample_count, block_align;
		int32_t sample_rate;
	
		while (!we_are_done) {
	
			struct {
				char chunk_type[4];
				int32_t chunk_length;
			} chunk_header;
	
			int32_t chunk_length, chunk_leftover;
	
			items_read = fread(&chunk_header, 8, 1, stdin);
	
	   	if (items_read != 1) {
				if (feof(stdin)) {
					fprintf(stderr, "Unexpected EOF\n");
					exit(1);
				} else {
		      	fprintf(stderr, "Header unreadable\n");
			      exit(1);
			   }
			}
	
			/* Now we have the chunk header with type and length */
	
			if (machine_endianness != file_endianness) {
				chunk_length = swap_int32(chunk_header.chunk_length);
			} else {
				chunk_length = chunk_header.chunk_length;
			}
	
			chunk_leftover = chunk_length;
	
			if (strncmp(chunk_header.chunk_type, "fmt ", 4) == 0) {
				if (debug_flag) fprintf(stderr, "Found fmt chunk with length %d\n", chunk_length);
	
				/* these are all little-endian */
				struct {
					int16_t audio_format;
					int16_t channel_count;
					int32_t sample_rate;
					int32_t byte_rate;
					int16_t block_align;
					int16_t bits_per_sample;
					int16_t extended_size;
					int16_t valid_bits_per_sample;
					int32_t channel_mask;
					int16_t ext_audio_format;
					char rest_of_subformat[14];
				} fmt_chunk;
	
				/* initialize this to 0, chunk may be shorter */
				fmt_chunk.extended_size = 0;

				if (chunk_length > sizeof(fmt_chunk)) {
					fprintf(stderr, "fmt chunk length is too long, got %d but expected %lu\n", chunk_length, sizeof(fmt_chunk));
					exit(1);
				}
	
				items_read = fread(&fmt_chunk, chunk_length, 1, stdin);
	
				if (items_read != 1) {
					fprintf(stderr, "Unreadable fmt chunk\n");
					exit(1);
				}
	
				int16_t audio_format = (machine_endianness != file_endianness ?  swap_int16(fmt_chunk.audio_format) : fmt_chunk.audio_format);
				int16_t extended_size = (machine_endianness != file_endianness ?  swap_int16(fmt_chunk.extended_size) : fmt_chunk.extended_size);

				if (audio_format == -2) {

					if (extended_size == 0) {
						fprintf(stderr, "Expecting extended fmt chunk, but it's missing\n");
						exit(1);
					}

					audio_format = (machine_endianness != file_endianness ?  swap_int16(fmt_chunk.ext_audio_format) : fmt_chunk.ext_audio_format);
				}
	
				if (audio_format != 1) {
					fprintf(stderr, "Unusable wave format tag: %d\n", audio_format);
					exit(1);
				}
	
				channel_count = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.channel_count) : fmt_chunk.channel_count);
				block_align = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.block_align) : fmt_chunk.block_align);
				bits_per_sample = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.bits_per_sample) : fmt_chunk.bits_per_sample);
	
				sample_rate = (machine_endianness != file_endianness ? swap_int32(fmt_chunk.sample_rate) : fmt_chunk.sample_rate);

				/* this is a WAV file, so we set the signed/unsigned flag here */
				file_signing = (bits_per_sample == 8 ? UNSIGNED : SIGNED);
	
				if (debug_flag) {
					fprintf(stderr, "Audio Format: %d, Channel Count: %d, Block Align: %d, Bits Per Sample: %d, Sample Rate: %d, Extended fmt Size: %d, Data Format: ",
						audio_format, channel_count, block_align, bits_per_sample, sample_rate, extended_size);
					fprintf(stderr, "%s\n", (file_signing == SIGNED ? "Signed" : "Unsigned"));
				}

				chunk_leftover = 0;
	
			} else if (strncmp(chunk_header.chunk_type, "data", 4) == 0) {
				if (debug_flag) fprintf(stderr, "Found data chunk with length %d\n", chunk_length);
				we_are_done = 1;
	
				sample_count = chunk_length / (bits_per_sample * channel_count / 8);
	
				if (debug_flag) fprintf(stderr, "Sample Count: %d\n", sample_count);
	
				calculate_waveform(stdin, sample_count, channel_count, bits_per_sample, algorithm, machine_endianness, data_endianness, file_signing);
	
			} else if (strncmp(chunk_header.chunk_type, "COMM", 4) == 0) {
				if (debug_flag) fprintf(stderr, "Found COMM chunk with length %d\n", chunk_length);
	
				/* these are all big-endian */
				struct {
					int16_t channel_count;
					int32_t sample_count;
					int16_t bits_per_sample;
					long double sample_rate;
					char padding[30]; /* just in case */
				} comm_chunk;
	
				if (chunk_length > sizeof(comm_chunk)) {
					fprintf(stderr, "comm chunk length is too long, got %d but expected %lu\n", chunk_length, sizeof(comm_chunk));
					exit(1);
				}
	
				items_read = fread(&comm_chunk, chunk_length, 1, stdin);
	
				if (items_read != 1) {
					fprintf(stderr, "Unreadable comm chunk\n");
					exit(1);
				}
	
				channel_count = (machine_endianness != file_endianness ? swap_int16(comm_chunk.channel_count) : comm_chunk.channel_count);
				bits_per_sample = (machine_endianness != file_endianness ? swap_int16(comm_chunk.bits_per_sample) : comm_chunk.bits_per_sample);
				sample_count = (machine_endianness != file_endianness ? swap_uint32(comm_chunk.sample_count) : comm_chunk.sample_count);
	
				if (debug_flag) {
					fprintf(stderr, "Channel Count: %d, Bits Per Sample: %d, Sample Count: %d, Data Format: ",
						channel_count, bits_per_sample, sample_count);
					fprintf(stderr, "%s\n", (file_signing == SIGNED ? "Signed" : "Unsigned"));
				}

				chunk_leftover = 0;
	
			} else if (strncmp(chunk_header.chunk_type, "SSND", 4) == 0) {
				if (debug_flag) fprintf(stderr, "Found SSND chunk with length %d\n", chunk_length);
				we_are_done = 1;
	
				/* The first 8 bytes are offset and blocksize */
				struct {
					uint32_t ssnd_offset, ssnd_blocksize;
				} comm_subheader;
	
				uint32_t ssnd_offset, ssnd_blocksize;
	
				items_read = fread(&comm_subheader, sizeof(comm_subheader), 1, stdin);
	
				if (items_read != 1) {
					fprintf(stderr, "Unreadable comm chunk\n");
					exit(1);
				}
	
				ssnd_offset = (machine_endianness != file_endianness ? swap_uint32(comm_subheader.ssnd_offset) : comm_subheader.ssnd_offset);
				ssnd_blocksize = (machine_endianness != file_endianness ? swap_uint32(comm_subheader.ssnd_blocksize) : comm_subheader.ssnd_blocksize);
	
				if (debug_flag) {
					fprintf(stderr, "ssnd offset: %d, ssnd blocksize: %d\n", ssnd_offset, ssnd_blocksize);
				}
	
				/* if there's an offset we'll skip that many bytes */
				if (ssnd_offset > 0) {
					fseek(stdin, ssnd_offset, SEEK_CUR);
				}
	
				calculate_waveform(stdin, sample_count, channel_count, bits_per_sample, algorithm, machine_endianness, data_endianness, file_signing);
	
				chunk_leftover = 0;
			}
	
			/* There will be a padding byte to skip if the chunk has an
			 * odd length */
			if ((chunk_length & 1) == 1) {
				chunk_leftover++;
			}

			if (fseek(stdin, chunk_leftover, SEEK_CUR) != 0) {
				if (feof(stdin)) {
					fprintf(stderr, "Unexpected EOF\n");
					exit(1);
				} else {
		      	fprintf(stderr, "Chunk unreadable\n");
			      exit(1);
			   }
			}
	
		}
	}
}
