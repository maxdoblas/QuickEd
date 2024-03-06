/*
 *                             The MIT License
 *
 * This file is part of QuickEdit library.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "align_benchmark_params.h"

/*
 * Default parameters
 */
align_bench_params_t parameters = {
  // Algorithm
  .algorithm = alignment_edit_quicked,
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
  // Other algorithms parameters
  .bandwidth = 15,
  .window_size = 9,
  .overlap_size = 1,
  .hew_percentage = 15,
  .hew_threshold = 40,
  .force_scalar = false,
  .only_score = false,
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 100000,
  .verbose = 0,
};
/*
 * Menu
 */
void usage(void) {
  fprintf(stderr,
      "USE: ./align_benchmark -a ALGORITHM -i PATH                             \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a ALGORITHM                                       \n"
      "            [Quicked Edit (Levenshtein)]                                \n"
      "              quicked                                                   \n"
      "              edit-banded                                               \n"
      "              edit-windowed                                             \n"
      "            [Other Edit (Levenshtein)]                                  \n"
      "              edlib                                                     \n"
      "              edit-dp                                                   \n"
      "              edit-dp-banded                                            \n"
      "              edit-bpm                                                  \n"
      "        [Input & Output]                                                \n"
      "          --input|i PATH                                                \n"
      "          --output|o PATH                                               \n"
      "          --output-full PATH                                            \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth INT                                               \n"
      "          --window-size INT                                             \n"
      "          --overlap-size INT                                            \n"
      "          --hew-threshold INT                                           \n"
      "          --hew-prercentage INT                                         \n"
      "          --only-score                                                  \n"
      "          --force-scalar                                                \n"
      "        [Misc]                                                          \n"
      "          --check|c 'display'|'correct'|'score'|'alignment'             \n"
      "          --check-bandwidth INT                                         \n"
      "        [System]                                                        \n"
      "          --num-threads|t INT                                           \n"
      "          --batch-size INT                                              \n"
      "          --progress|P INT                                              \n"
      "          --verbose|v INT                                               \n"
      "          --quiet|q                                                     \n"
      "          --help|h                                                      \n");
}
/*
 * Parsing arguments
 */
void parse_arguments(
    int argc,
    char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
    { "window-size", required_argument, 0, 2001 },
    { "overlap-size", required_argument, 0, 2002 },
    { "hew-threshold", required_argument, 0, 2003 },
    { "hew-percentage", required_argument, 0, 2004 },
    { "force-scalar", no_argument, 0, 2005 },
    { "only-score", no_argument, 0, 2006 },
    /* Misc */
    { "check", required_argument, 0, 'c' },
    { "check-bandwidth", required_argument, 0, 3002 },
    /* System */
    { "num-threads", required_argument, 0, 't' },
    { "batch-size", required_argument, 0, 4000 },
    { "progress", required_argument, 0, 'P' },
    { "verbose", required_argument, 0, 4001 },
    { "verbose1", no_argument, 0, 'v' },
    { "quiet", no_argument, 0, 'q' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:vqt:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Algorithm
     */
    case 'a': {
      // Edit
      if (strcmp(optarg,"edit-bpm")==0) {
        parameters.algorithm = alignment_edit_bpm;
      } else if (strcmp(optarg,"edit-banded")==0) {
        parameters.algorithm = alignment_edit_banded;
      } else if (strcmp(optarg,"edit-banded-hirschberg")==0) {
        parameters.algorithm = alignment_edit_banded_hirschberg;
      } else if (strcmp(optarg,"quicked")==0) {
        parameters.algorithm = alignment_edit_quicked;
      } else if (strcmp(optarg,"edit-windowed")==0) {
        parameters.algorithm = alignment_edit_windowed;
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
      }  else if (strcmp(optarg,"edlib")==0) {
        parameters.algorithm = alignment_edlib;
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    /*
     * Input/Output
     */
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    case 800: // --output-full
      parameters.output_filename = optarg;
      parameters.output_full = true;
      break;
    /*
     * Other alignment parameters
     */
    case 2000: // --bandwidth
      parameters.bandwidth = atoi(optarg);
      break;
    case 2001: // --window-size
      parameters.window_size = atoi(optarg);
      break;
    case 2002: // --overlap-size
      parameters.overlap_size = atoi(optarg);
      break;
    case 2003: // --hew-threshold
      parameters.hew_threshold = atoi(optarg);
      break;
    case 2004: // --hew-percentage
      parameters.hew_percentage = atoi(optarg);
      break;
    case 2005: // --force-scalar
      parameters.force_scalar = true;
      break;
    case 2006: // --only-score
      parameters.only_score = true;
      break;
    /*
     * Misc
     */
    case 'c':
      if (optarg ==  NULL) { // default = correct
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"display")==0) {
        parameters.check_display = true;
      } else if (strcasecmp(optarg,"correct")==0) {
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"score")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"alignment")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = true;
      } else {
        fprintf(stderr,"Option '--check' must be in {'correct','score','alignment'}\n");
        exit(1);
      }
      break;
    case 3002: // --check-bandwidth
      parameters.check_bandwidth = atoi(optarg);
      break;
    /*
     * System
     */
    case 't': // --num-threads
      parameters.num_threads = atoi(optarg);
      break;
    case 4000: // --batch-size
      parameters.batch_size = atoi(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      parameters.verbose = 1;
      break;
    case 4001: // --verbose (long option)
      parameters.verbose = atoi(optarg);
      if (parameters.verbose < 0 || parameters.verbose > 4) {
        fprintf(stderr,"Option '--verbose' must be in {0,1,2,3,4}\n");
        exit(1);
      }
      break;
    case 'q':
      parameters.verbose = -1;
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // Checks general
  if (parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  // Check 'bandwidth' parameter
  switch (parameters.algorithm) {
    case alignment_edit_dp_banded:
    case alignment_edit_banded:
    case alignment_edit_banded_hirschberg:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      } else if (parameters.bandwidth < 1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be > 0\n");
        exit(1);
      }
      break;
    case alignment_edit_windowed:
      if (parameters.window_size == -1) {
        fprintf(stderr,"Parameter 'window-size' has to be provided for banded algorithms\n");
        exit(1);
      } else if (parameters.window_size < 1) {
        fprintf(stderr,"Parameter 'window-size' has to be > 0\n");
        exit(1);
      }

      if (parameters.overlap_size == -1) {
        fprintf(stderr,"Parameter 'overlap-size' has to be provided for banded algorithms\n");
        exit(1);
      } else if (parameters.overlap_size > parameters.window_size - 1 || parameters.overlap_size < 0) {
        fprintf(stderr,"Parameter 'overlap-size' has to be: 0 <= overlap-size < window-size \n");
        exit(1);
      }
      break;
    case alignment_edit_dp:
    case alignment_edit_bpm:
    case alignment_edlib:
    case alignment_edit_quicked:
    default:
      break;
  }
}
