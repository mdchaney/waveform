/*
	waveform

	By Michael Chaney
	Copyright 2014 Michael Chaney Consulting Corporation - All rights reserved

	Released under terms of the MIT license

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
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

static int verbose_flag;
static unsigned long points;
static int use_mean;
static int use_rms;

typedef enum { BIG, LITTLE, MIXED } Endianness_t;
typedef enum { RIFF_FILE, FORM_FILE } FileFormat_t;
typedef enum { WAVE_FILE, AIFF_FILE, AIFC_FILE } AudioFormat_t;

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

int main(int argc, char **argv) {
	int c;

	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 },
			{ "mean", no_argument, &use_mean, 1 },
			{ "rms", no_argument, &use_rms, 1 },
			{ "points", required_argument, 0, 'p' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		c = getopt_long(argc, argv, "p", long_options, &option_index);

		if (c == -1) break;

		switch(c) {
			case 0:
				if (long_options[option_index].flag != 0) break;
				break;

			case 'p':
				/* parse optarg to get # */
				if (sscanf(optarg, "%lu", &points) != 1) {
					printf("Expecting number for 'points', got `%s'\n", optarg);
					exit(1);
				}
				break;

			case '?': /* error condition */
				break;

			default:
				abort();
		}
	}

	if (!use_mean && !use_rms) use_rms = 1;

	if (verbose_flag) {
		printf("verbose!\n");
	}

	if (points) {
		printf("points: %lu\n", points);
	}

	if (use_mean) {
		printf("using mean\n");
	} else if (use_rms) {
		printf("using rms\n");
	} else {
		printf("confused about algorithm\n");
	}

	char* filename;

	fprintf(stderr, "optind %d, argc %d\n", optind, argc);

	if (optind-1 < argc) {
		filename = argv[optind];
	}

	if (filename) {
		printf("Filename: %s\n", filename);
	}

	if (!filename) {
		printf("Reading stdin\n");
	}

	if (filename) {
		if (!freopen(filename, "r", stdin)) {
			printf("Cannot open `%s' for reading\n", filename);
			exit(1);
		}
	}

	Endianness_t machine_endianness, file_endianness, data_endianness;

	uint32_t i=0x01234567;

	switch (*((uint8_t*)(&i))) {
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

	printf("Your architecture is ");

	if (machine_endianness == BIG) {
		printf("big endian\n");
	} else if (machine_endianness == LITTLE) {
		printf("little endian\n");
	} else if (machine_endianness == MIXED) {
		printf("mixed endian\n");
	}

	size_t items_read;

	/* At this point, options are parsed and the file is opened.  We will
	 * read the first 8 bytes and see if it's a RIFF header. */

	typedef struct RIFF_HEADER {
		char RIFF[4];
		int32_t file_size_sans_riff_header;
	} riff_header_struct;

	riff_header_struct riff_header;

	items_read = fread(&riff_header, sizeof(riff_header_struct), 1, stdin);

	if (items_read != 1) {
		printf("Header unreadable\n");
		exit(1);
	}

	if (strncmp(riff_header.RIFF, "RIFF", 4) != 0 && strncmp(riff_header.RIFF, "FORM", 4) != 0) {
		printf("No RIFF/FORM header\n");
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
		printf("FORM file format\n");
		file_endianness = BIG;
	} else if (file_format == RIFF_FILE) {
		printf("RIFF file format\n");
		file_endianness = LITTLE;
	} else {
		printf("Confusion in file format\n");
		exit(1);
	}

	printf("Your file is ");

	if (file_endianness == BIG) {
		printf("big endian\n");
	} else if (file_endianness == LITTLE) {
		printf("little endian\n");
	}

	if (file_endianness != machine_endianness) {
		printf("File size (sans RIFF/FORM header): %d\n", swap_int32(riff_header.file_size_sans_riff_header));
	} else {
		printf("File size (sans RIFF/FORM header): %d\n", riff_header.file_size_sans_riff_header);
	}

	/* Now, we need to see what kind of file it is.  It's a RIFF if it
	 * got this far, this should be WAVE, AIFF or AIFC. */

	char audio_type[4];

	items_read = fread(audio_type, 4, 1, stdin);

	AudioFormat_t audio_format;

	printf("File type: %c%c%c%c\n", audio_type[0], audio_type[1], audio_type[2], audio_type[3]);

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
			fprintf(stderr, "Found fmt chunk with length %d\n", chunk_length);

			/* these are all little-endian */
			struct {
				int16_t audio_format;
				int16_t channel_count;
				int32_t sample_rate;
				int32_t byte_rate;
				int16_t block_align;
				int16_t bits_per_sample;
				char    padding[20]; /* sometimes the chunk is longer */
			} fmt_chunk;

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

			if (audio_format != 1) {
				fprintf(stderr, "Unusable wave format tag: %d\n", audio_format);
				exit(1);
			}

			channel_count = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.channel_count) : fmt_chunk.channel_count);
			block_align = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.block_align) : fmt_chunk.block_align);
			bits_per_sample = (machine_endianness != file_endianness ? swap_int16(fmt_chunk.bits_per_sample) : fmt_chunk.bits_per_sample);

			sample_rate = (machine_endianness != file_endianness ? swap_int32(fmt_chunk.sample_rate) : fmt_chunk.sample_rate);

			fprintf(stderr, "Audio Format: %d, Channel Count: %d, Block Align: %d, Bits Per Sample: %d, Sample Rate: %d\n",
				audio_format, channel_count, block_align, bits_per_sample, sample_rate);

			chunk_leftover = 0;

		} else if (strncmp(chunk_header.chunk_type, "data", 4) == 0) {
			fprintf(stderr, "Found data chunk with length %d\n", chunk_length);
			we_are_done = 1;

			sample_count = chunk_length / (bits_per_sample * channel_count / 8);

			fprintf(stderr, "Sample Count: %d\n", sample_count);

			/* Here's where we implement Bresenham's algorithm for line
			 * drawing.  Let's assume a line in the first octant (in polar
			 * coordinates the angle is 0-45 degrees) where "x" is the
			 * actual sample_count and "y" is the "points" variable.  We
			 * assume then that "sample_count > points".
			 *
			 * In standard Bresenham talk that means that "sample_count" is
			 * "dx" and "points" is "dy".  Using integer division the
			 * number of samples per point will vary between two numbers,
			 * say "j" and "j+1".  When to use "j" and when to use "j+1" is
			 * the trick.
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
			 * use "j" samples and when we use "j+1".  As the "slope"
			 * exceeds 45 degrees (at "points / 2") we have to use "j+1"
			 * more often than "j".
			 *
			 * Take the first case, where samples is between 9000 and 9499.
			 * Generalized (with integer arithmetic):
			 *
			 * j = samples / points
			 * leftover = samples - (j * points)
			 *
			 * At this point with the example numbers above "samples" is
			 * "1000" * and "leftover" is the last column.  For our example
			 * we examine numbers where "leftover / points < 0.5" for
			 * simplicity.  In that case we need to stick between 0 and 499 
			 * "j+1" sample groups in with the other samples.  So every so
			 * often we have to grab "j+1" instead of "j".
			 *
			 */

			long j = sample_count / points;
			long leftover = sample_count - (j * points);

			fprintf(stderr, "j: %ld, leftover: %ld\n", j, leftover);

			long samples_left = sample_count;
			int sample_group_sizes[10000];

			int i = 0, jump_at, jump_counter = 0;
			int under_45=1;

			if (leftover > points / 2) {
				under_45 = 0;
				jump_at = points / (points - leftover);

				while (samples_left > 0) {
					if (jump_counter < jump_at) {
						samples_left -= j+1;
						sample_group_sizes[i] = j+1;
					} else {
						samples_left -= j;
						sample_group_sizes[i] = j;
						jump_counter = 0;
					}
					i++;
					jump_counter++;
				}

			} else {
				under_45 = 1;
				jump_at = points / (leftover + 1);

				while (samples_left > 0) {
					if (jump_counter < jump_at) {
						samples_left -= j;
						sample_group_sizes[i] = j;
					} else {
						samples_left -= j+1;
						sample_group_sizes[i] = j+1;
						jump_counter = 0;
					}
					i++;
					jump_counter++;
				}

			}

			samples_left = sample_count;

			for (i=0 ; i<points; i++) {
				fprintf(stderr, "%5d:   %10d   %10d\n", i, sample_group_sizes[i], samples_left);
				samples_left -= sample_group_sizes[i];
			}

			fprintf(stderr, "%5d:                %10d\n", i, samples_left);

		} else if (strncmp(chunk_header.chunk_type, "COMM", 4) == 0) {
			fprintf(stderr, "Found COMM chunk with length %d\n", chunk_length);
		} else if (strncmp(chunk_header.chunk_type, "SSND", 4) == 0) {
			fprintf(stderr, "Found SSND chunk with length %d\n", chunk_length);
			we_are_done = 1;
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
