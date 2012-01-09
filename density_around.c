/*
 density_around.c
 project: DensityAround
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 06-Jan-2012

 see README

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "string_buffer.h"
#include "utility_lib.h"

typedef unsigned long bin_t;

// For this run
char* cmd;

// Command line arguments
unsigned long bin_size, num_of_bins;
char *objects_file_path, *events_file_path;

// Optional cmdline args
char *bin_denom_path = NULL;
FILE* bin_count_out = NULL;
long region_start = 0, region_end = 0;

// Events data
unsigned long num_of_events = 0;
unsigned long events_arr_capacity = 500000;
long *event_starts, *event_ends;
unsigned long sum_event_lengths = 0;

double mean_event_length;

// Bin results
bin_t num_events_overlapping_objects = 0;
bin_t *bins_left, *bins_right;

// Denominator results
double object_overlap_bin = 0;
double *bins_left_denom = NULL, *bins_right_denom = NULL;
unsigned long bins_left_denom_all = 0, bins_right_denom_all = 0;

void print_usage(char* err_msg)
{
  if(err_msg != NULL)
  {
    fprintf(stderr, "Error: %s\n", err_msg);
  }

  fprintf(stderr, "usage: %s [OPTIONS] <bin_size> <num_of_bins> <events.csv> "
                  "<objects.csv> <out.csv>\n", cmd);
  fprintf(stderr, "  Histogram of density of events around objects\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  - objects and events must be sorted csv files "
                  "each line reading: 'start,end'\n");
  fprintf(stderr, "  - bin_size is the width of the bins\n");
  fprintf(stderr, "  - num_of_bins is the number of bins either side\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  OPTIONS:\n");
  fprintf(stderr, "    --density <start> <end> <denominators.csv>\n");
  fprintf(stderr, "    If specified will give counts of bin occurance, required\n");
  fprintf(stderr, "    to calculate actual event density around objects (e.g.\n");
  fprintf(stderr, "    '--density 1 2000000 chr1.csv' for a chromosome)\n");
  fprintf(stderr, "\n");
  
  fprintf(stderr, "  turner.isaac@gmail.com  06 Jan 2012\n");

  exit(EXIT_FAILURE);
}

// Update bin denominators
inline void update_bin_denoms(unsigned long *all_bins_denom, double *bins_denom,
                              const long remaining)
{
  unsigned long last_bin_index = remaining / bin_size;

  if(last_bin_index >= num_of_bins)
  {
    (*all_bins_denom)++;
  }
  else
  {
    unsigned long i;
    for(i = 0; i < last_bin_index; i++)
    {
      bins_denom[i]++;
    }

    bins_denom[last_bin_index] += (double)(remaining % bin_size) / bin_size;
  }
}

void load_all_events(gzFile* events_file, char reading_from_stdin)
{
  STRING_BUFFER* events_line = string_buff_init(200);

  t_buf_pos read_length;
  unsigned long line_num;
  char seen_header = 0;
  char prev_line_emtpy = 0;

  for(line_num = 1;
      (read_length = string_buff_reset_readline(events_line, events_file)) > 0;
      line_num++)
  {
    string_buff_chomp(events_line);

    // Look for two empty lines in a row to mark end of events input
    if(reading_from_stdin)
    {
      if(strlen(events_line->buff) == 0)
      {
        if(prev_line_emtpy == 1)
        {
          break;
        }

        prev_line_emtpy = 1;
      }
      else
      {
        prev_line_emtpy = 0;
      }
    }

    //printf("reading events %lu '%s'\n", read_length, events_line->buff);

    // Check if comment line
    if(events_line->buff[0] != '#' &&
       !string_is_all_whitespace(events_line->buff))
    {
      char* separator = strchr(events_line->buff, ',');
      
      if(separator == NULL)
      {
        separator = strchr(events_line->buff, '\t');
      }

      if(separator == NULL)
      {
        fprintf(stderr, "Error (file: %s line: %lu): no separator found "
                        "(tab or comma only)\n",
                objects_file_path, line_num);
        print_usage(NULL);
      }

      *separator = '\0';

      char* start_num_str = trim(events_line->buff);
      char* end_num_str = trim(separator+1);

      //printf("%s:%lu - '%s' , '%s'\n", events_file_path, line_num,
      //       start_num_str, end_num_str);

      long event_start, event_end;

      char successs = (parse_entire_long(start_num_str, &event_start) &&
                       parse_entire_long(end_num_str, &event_end) &&
                       event_end >= event_start);

      if(!successs)
      {
        if(num_of_events > 0 || seen_header)
        {
          // Error
          fprintf(stderr, "Error on line: %lu file: %s\n",
                  line_num, events_file_path);
          print_usage(NULL);
        }

        // Everyone gets one
        seen_header = 1;
      }
      else
      {
        sum_event_lengths += event_end - event_start;

        if(num_of_events + 1 == events_arr_capacity)
        {
          // Need to expand array
          events_arr_capacity *= 2;
          event_starts = realloc(event_starts, events_arr_capacity * sizeof(long));
          event_ends = realloc(event_ends, events_arr_capacity * sizeof(long));
      
          if(event_starts == NULL || event_ends == NULL)
          {
            print_usage("Not enough memory to read in events!");
          }
        }

        event_starts[num_of_events] = event_start;
        event_ends[num_of_events] = event_end;
        num_of_events++;
      }
    }
  }

  string_buff_free(events_line);
}

void run_through_objects(gzFile* objects_file)
{
  STRING_BUFFER* objects_line = string_buff_init(200);

  t_buf_pos read_length;
  unsigned long line_num;
  char seen_header = 0;

  for(line_num = 1;
      (read_length = string_buff_reset_readline(objects_line, objects_file)) > 0;
      line_num++)
  {
    string_buff_chomp(objects_line);
    //printf("reading objects %lu '%s'\n", read_length, objects_line->buff);

    // Check if comment line
    if(objects_line->len > 0 &&
       objects_line->buff[0] != '#' &&
       !string_is_all_whitespace(objects_line->buff))
    {
      char* separator = strchr(objects_line->buff, ',');
      
      if(separator == NULL)
      {
        separator = strchr(objects_line->buff, '\t');
      }

      if(separator == NULL)
      {
        fprintf(stderr, "Error (file: %s line: %lu): no separator found "
                        "(tab or comma only)\n",
                objects_file_path, line_num);
        print_usage(NULL);
      }

      *separator = '\0';

      char* start_num_str = trim(objects_line->buff);
      char* end_num_str = trim(separator+1);

      long obj_start, obj_end;

      char successs = (parse_entire_long(start_num_str, &obj_start) &&
                       parse_entire_long(end_num_str, &obj_end) &&
                       obj_start > 0 && obj_start <= obj_end);

      if(!successs && seen_header)
      {
        // Error
        fprintf(stderr, "Error (file: %s line: %lu): not valid numbers\n",
                objects_file_path, line_num);
        print_usage(NULL);
      }

      seen_header = 1;

      if(successs)
      {
        // Read in successfully
      
        if(bin_count_out != NULL)
        {
          // update bins with this object

          if(obj_start < region_start)
          {
            fprintf(stderr, "Warning: object starts before region (%li < %li)\n",
                    obj_start, region_start);
          }
          else if(obj_end > region_end)
          {
            fprintf(stderr, "Warning: object ends after region (%li > %li)\n",
                    obj_end, region_end);
          }
          else
          {
            object_overlap_bin += obj_end - obj_start + 2 * mean_event_length;

            long left_remaining = obj_start - region_start;

            update_bin_denoms(&bins_left_denom_all, bins_left_denom,
                              left_remaining);

            long right_remaining = region_end - obj_end;

            update_bin_denoms(&bins_right_denom_all, bins_right_denom,
                              right_remaining);
          }
        } // done updating bin denominators

        // Loop through events, updating bins using position start, end
        unsigned long events_i;
      
        for(events_i = 0; events_i < num_of_events; events_i++)
        {
          long event_start = event_starts[events_i];
          long event_end = event_ends[events_i];

          if((event_start >= obj_start && event_start <= obj_end) ||
             (event_end >= obj_start && event_end <= obj_end))
          {
            num_events_overlapping_objects++;
          }
          else if(event_end < obj_start)
          {
            // Event to the left of the object
            long dist = event_start - obj_end;
            bin_t bin = dist % bin_size;

            if(bin < num_of_bins)
            {
              bins_left[bin]++;
            }
          }
          else
          {
            // Event to the right of the object
            long dist = obj_start - event_end;
            bin_t bin = dist % bin_size;

            if(bin < num_of_bins)
            {
              bins_right[bin]++;
            }
          }
        } // loop over events
      } // if(success)
    } // if(!comment && !whitespace)
  } // loop over objects

  string_buff_free(objects_line);
}

int main(int argc, char* argv[])
{
  cmd = argv[0];

  #ifdef DEBUG
  printf("DEBUG: on\n");
  #endif

  if(argc != 6 && argc != 10)
  {
    print_usage(NULL);
  }

  // Parse command line args
  int argi = 1;

  if(argc == 10)
  {
    if(strcasecmp(argv[1], "--density") != 0)
    {
      fprintf(stderr, "Error: unknown arg (expected --length) '%s'\n", argv[1]);
      print_usage(NULL);
    }
    else if(!parse_entire_long(argv[2], &region_start))
    {
      fprintf(stderr, "Error: Invalid region start arg '%s'\n", argv[2]);
      print_usage(NULL);
    }
    else if(!parse_entire_long(argv[3], &region_end))
    {
      fprintf(stderr, "Error: Invalid region end arg '%s'\n", argv[3]);
      print_usage(NULL);
    }
    else if(region_start >= region_end)
    {
      print_usage("region start must be less than the region end\n");
    }

    bin_denom_path = argv[4];

    if(strcmp(bin_denom_path,"-") == 0)
    {
      bin_count_out = stdout;
    }
    else
    {
      bin_count_out = fopen(bin_denom_path, "w");
    }

    if(bin_count_out == NULL)
    {
      fprintf(stderr, "Error: Cannot open output file '%s'\n", argv[4]);
      print_usage(NULL);
    }

    // start reading arguments from argv[5]
    argi = 5;
  }
  else
  {
    // Start reading arguments from argv[1]
    argi = 1;
  }
  
  char* bin_size_arg = argv[argi++];
  char* num_of_bins_arg = argv[argi++];
  objects_file_path = argv[argi++];
  events_file_path = argv[argi++];
  char* output_path = argv[argi++];

  // Parse bin size, max bin numbers
  if(!parse_entire_ulong(bin_size_arg, &bin_size) || bin_size <= 0)
  {
    print_usage("Invalid <bin_size> argument");
  }

  if(!parse_entire_ulong(num_of_bins_arg, &num_of_bins) || num_of_bins <= 0)
  {
    print_usage("Invalid <num_of_bins> argument");
  }

  // Open input files
  gzFile *stdin_file = NULL, *events_file = NULL, *objects_file = NULL;

  if(strcmp(events_file_path, "-") == 0 && strcmp(objects_file_path ,"-") == 0)
  {
    // Open stdin
    stdin_file = gzdopen(fileno(stdin), "r");
  }
  else
  {
    events_file = gzopen(events_file_path, "r");
    objects_file = gzopen(objects_file_path, "r");
  
    if(events_file == NULL)
    {
      fprintf(stderr, "Error: Cannot open events file '%s'\n", events_file_path);
      print_usage(NULL);
    }

    if(objects_file == NULL)
    {
      fprintf(stderr, "Error: Cannot open objects file '%s'\n", objects_file_path);
      print_usage(NULL);
    }
  }

  // Open output file
  FILE* out;
  
  if(strcmp(output_path, "-") == 0)
  {
    out = stdout;
  }
  else
  {
    out = fopen(output_path, "w");
  }


  //
  // Load all events
  //
  
  // Initialise array to hold events
  num_of_events = 0;
  events_arr_capacity = 500000;
  event_starts = (long*) malloc(events_arr_capacity * sizeof(long));
  event_ends = (long*) malloc(events_arr_capacity * sizeof(long));

  sum_event_lengths = 0;

  if(event_starts == NULL || event_ends == NULL)
  {
    print_usage("Not enough memory to read in events!");
  }

  if(stdin_file != NULL)
  {
    load_all_events(stdin_file, 1);
  }
  else
  {
    load_all_events(events_file, 0);
  }

  mean_event_length = (num_of_events == 0) ? 0 : sum_event_lengths / num_of_events;

  if(out != stdout || 1)
  {
    printf("Loaded %lu events\n", num_of_events);
    printf("Mean event length %f\n", mean_event_length);
  }

  // Create bins
  num_events_overlapping_objects = 0;
  bins_left = (bin_t*) malloc(num_of_bins * sizeof(bin_t));
  bins_right = (bin_t*) malloc(num_of_bins * sizeof(bin_t));

  if(bins_left == NULL || bins_right == NULL)
  {
    print_usage("Not enough memory to allocate bin arrays!");
  }

  // Initialise bins to zero
  unsigned long i;
  for (i = 0; i < num_of_bins; i++)
  {
    bins_left[i] = 0;
    bins_right[i] = 0;
  }

  if(bin_count_out != NULL)
  {
    bins_left_denom = (double*) malloc(num_of_bins * sizeof(double));
    bins_right_denom = (double*) malloc(num_of_bins * sizeof(double));

    for(i = 0; i < num_of_bins; i++)
    {
      bins_left_denom[i] = 0;
      bins_right_denom[i] = 0;
    }
  }

  //
  // Read objects one at a time 
  //
  
  run_through_objects(objects_file != NULL ? objects_file : stdin_file);

  if(out != stdout || 1)
  {
    printf("Saving bin counts to file: %s\n", output_path);
  }

  // Have to use < ULONG_MAX instead of >=0 since we're using unsigned longs
  for (i = num_of_bins-1; i < ULONG_MAX; i--)
  {
    fprintf(out, "%f,%lu\n", -(i+0.5)*bin_size, bins_left[i]);
  }

  fprintf(out, "0,%lu\n", num_events_overlapping_objects);

  for (i = 0; i < num_of_bins; i++)
  {
    fprintf(out, "%f,%lu\n", (i+0.5)*bin_size, bins_right[i]);
  }

  if(out != stdout)
  {
    fclose(out);
  }

  if(bin_count_out != NULL)
  {
    // Save bin denominators

    if(bin_count_out == stdout)
    {
      // Print separator
      printf("\n\n");
    }
    else
    {
      printf("Saving bin denominators to file: %s\n", bin_denom_path);
    }

    // Have to use < ULONG_MAX instead of >=0 since we're using unsigned longs
    for (i = num_of_bins-1; i < ULONG_MAX; i--)
    {
      fprintf(bin_count_out, "%f,%f\n", -(i+0.5)*bin_size,
              bins_left_denom[i] + bins_left_denom_all);
    }

    fprintf(bin_count_out, "0,%f\n", object_overlap_bin);

    for (i = 0; i < num_of_bins; i++)
    {
      fprintf(bin_count_out, "%f,%f\n", (i+0.5)*bin_size,
              bins_right_denom[i] + bins_right_denom_all);
    }
  }

  if(stdin_file != NULL)
  {
    gzclose(stdin_file);
  }
  else
  {
    gzclose(events_file);
    gzclose(objects_file);
  }

  return EXIT_SUCCESS;
}
