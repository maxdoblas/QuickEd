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
  .algorithm = alignment_edit_dp,
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
  // Alignment form
  .endsfree = false,
  .pattern_begin_free = 0.0,
  .text_begin_free = 0.0,
  .pattern_end_free = 0.0,
  .text_end_free = 0.0,
  // Other algorithms parameters
  .bandwidth = -1,
  .window_size = -1,
  .overlap_size = -1,
  .window_config = WINDOW_ALIGNED,
#ifdef EXTERNAL_BENCHMARKS
  /* ... */
#endif
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .plot = 0,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 100000,
  .verbose = 0,
};
/*
 * Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -a ALGORITHM -i PATH                             \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a ALGORITHM                                       \n"
      "            [Edit (Levenshtein)]                                        \n"
      "              edit-bpm                                                  \n"
      "              edit-bpm-banded                                           \n"
      "              edit-bpm-windowed                                         \n"
      "              edit-dp                                                   \n"
      "              edit-dp-banded                                            \n"
#ifdef EXTERNAL_BENCHMARKS
      "            [External/BitPal]                                           \n"
      "              bitpal-edit          (Edit)[score-only]                   \n"
      "            [External/Daligner]                                         \n"
      "              daligner             (Edit)                               \n"
      "            [External/Diffutils]                                        \n"
      "              diffutils            (Edit)                               \n"
      "            [External/Edlib]                                            \n"
      "              edlib                (Edit)                               \n"
      "            [External/LV89]                                             \n"
      "              lv89                 (Edit)[score-only]                   \n"
      "            [External/Parasail]                                         \n"
      "              parasail-nw-stripped (Gap-affine)                         \n"
      "              parasail-nw-scan     (Gap-affine)                         \n"
      "              parasail-nw-diag     (Gap-affine)                         \n"
      "              parasail-nw-banded   (Gap-affine)[score-only]             \n"
      "            [External/Scrooge]                                          \n"
      "              scrooge              (Edit)                               \n"
      "            [External/SeqAn]                                            \n"
      "              seqan-edit           (Edit)                               \n"
      "              seqan-edit-bpm       (Edit)[score-only]                   \n"
#endif
      "        [Input & Output]                                                \n"
      "          --input|i PATH                                                \n"
      "          --output|o PATH                                               \n"
      "          --output-full PATH                                            \n"
      "        [Penalties & Span]                                              \n"
      "          --ends-free P0,Pf,T0,Tf                                       \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth INT                                               \n"
#ifdef EXTERNAL_BENCHMARKS
      /* ... */
#endif
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-bandwidth INT                                         \n"
    //"          --plot                                                        \n"
      "        [System]                                                        \n"
    //"          --num-threads|t INT                                           \n"
    //"          --batch-size INT                                              \n"
    //"          --progress|P INT                                              \n"
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
    /* Penalties */
    { "ends-free", required_argument, 0, 901 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
    { "window-size", required_argument, 0, 2001 },
    { "overlap-size", required_argument, 0, 2002 },
    { "window-config", required_argument, 0, 2003 },
#ifdef EXTERNAL_BENCHMARKS
    /* ... */
#endif
    /* Misc */
    { "check", required_argument, 0, 'c' },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
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
      /* Test bench */
      if (strcmp(optarg,"test")==0) {
        parameters.algorithm = alignment_test;
      // Edit
      } else if (strcmp(optarg,"edit-bpm")==0) {
        parameters.algorithm = alignment_edit_bpm;
      } else if (strcmp(optarg,"edit-bpm-banded")==0) {
        parameters.algorithm = alignment_edit_bpm_banded;
      } else if (strcmp(optarg,"quicked")==0) {
        parameters.algorithm = alignment_edit_bpm_quicked;
      } else if (strcmp(optarg,"edit-bpm-windowed")==0) {
        parameters.algorithm = alignment_edit_bpm_windowed;
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
#ifdef EXTERNAL_BENCHMARKS
      /*
       * External Algorithm
       */
      // External (BitPal)
      } else if (strcmp(optarg,"bitpal-edit")==0) {
        parameters.algorithm = alignment_bitpal_edit;
      // External (Daligner)
      } else if (strcmp(optarg,"daligner")==0) {
        parameters.algorithm = alignment_daligner;
      // External (Diffutils)
      } else if (strcmp(optarg,"diffutils")==0) {
        parameters.algorithm = alignment_diffutils;
      // External (Edlib)
      } else if (strcmp(optarg,"edlib")==0) {
        parameters.algorithm = alignment_edlib;
      // External (LV89)
      } else if (strcmp(optarg,"lv89")==0) {
        parameters.algorithm = alignment_lv89;
      // External (Parasail)
      } else if (strcmp(optarg,"parasail-nw-stripped")==0) {
        parameters.algorithm = alignment_parasail_nw_stripped;
      } else if (strcmp(optarg,"parasail-nw-scan")==0) {
        parameters.algorithm = alignment_parasail_nw_scan;
      } else if (strcmp(optarg,"parasail-nw-diag")==0) {
        parameters.algorithm = alignment_parasail_nw_diag;
      } else if (strcmp(optarg,"parasail-nw-banded")==0) {
        parameters.algorithm = alignment_parasail_nw_banded;
      // External (Scrooge)
      } else if (strcmp(optarg,"scrooge")==0) {
        parameters.algorithm = alignment_scrooge;
      // External (SeqAn)
      } else if (strcmp(optarg,"seqan-edit")==0) {
        parameters.algorithm = alignment_seqan_edit;
      } else if (strcmp(optarg,"seqan-edit-bpm")==0) {
        parameters.algorithm = alignment_seqan_edit_bpm;
#endif
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
     * Penalties
     */
    case 901: { // --ends-free P0,Pf,T0,Tf
      parameters.endsfree = true;
      char* sentinel = strtok(optarg,",");
      parameters.pattern_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.pattern_end_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_begin_free = atof(sentinel);
      sentinel = strtok(NULL,",");
      parameters.text_end_free = atof(sentinel);
      break;
    }
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
    case 2003: // --overlap-size
      if (strcmp(optarg,"aligned")==0) {
        parameters.window_config = WINDOW_ALIGNED;
      } else if (strcmp(optarg,"unaligned")==0) {
        parameters.window_config = WINDOW_UNALIGNED;
      } else if (strcmp(optarg,"sse")==0) {
        parameters.window_config = WINDOW_SSE;
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
#ifdef EXTERNAL_BENCHMARKS
      /* ... */
#endif
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
    case 3003: // --plot
      parameters.plot = (optarg==NULL) ? 1 : atoi(optarg);
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
  if (parameters.algorithm!=alignment_test && parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  // Check 'bandwidth' parameter
  switch (parameters.algorithm) {
    case alignment_edit_bpm_banded:
    case alignment_edit_dp_banded:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      }
      if (parameters.window_size != -1) {
        fprintf(stderr,"Parameter 'window-size' has no effect with the selected algorithm\n");
        exit(1);
      }
      if (parameters.overlap_size != -1) {
        fprintf(stderr,"Parameter 'overlap-size' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
    case alignment_edit_bpm_windowed:
      if (parameters.bandwidth != -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }

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
    default:
      if (parameters.bandwidth != -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }
      if (parameters.window_size != -1) {
        fprintf(stderr,"Parameter 'window-size' has no effect with the selected algorithm\n");
        exit(1);
      }
      if (parameters.overlap_size != -1) {
        fprintf(stderr,"Parameter 'overlap-size' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
  }
  // Checks parallel
  if (parameters.num_threads > 1) {
    if (parameters.plot > 0) {
      fprintf(stderr,"Parameter 'plot' disabled for parallel executions\n");
      parameters.plot = 0;
    }
  }
}
