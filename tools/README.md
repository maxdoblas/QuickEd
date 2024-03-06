# TOOLS

## WHAT IS THIS?

Alongside the QuickEd library, we bundle tools that exploit the different features and serve as examples of how to use and exploit the library functions.
Also, these tools serve to test the library and benchmark its performance against other state-of-the-art tools and libraries for pairwise alignment.

* [Generate Dataset](#generate-dataset-tool)
* [Align Benchmark](#alignment-benchmark-tool)

## GENERATE DATASET TOOL

The *generate-dataset* tool allows the generating of synthetic random datasets for testing and benchmarking purposes.
This tool produces a simple output format (i.e., each pair of sequences in 2 lines) containing the pairs of sequences to be aligned using, for example, the [*align-benchmark*](#alignment-benchmark-tool) tool.
For example, the following command generates a dataset named 'sample.dataset.seq' of 5M pairs of 100 bases with an alignment error of 5% (i.e., 5 mismatches, insertions, or deletions per alignment).

```bash
./bin/generate_dataset -n 5000000 -l 100 -e 0.05 -o sample.dataset.seq
```

### Command-line Options (generate-dataset)

```bash
  --output|o <File>
    # Filename/Path to the output dataset.

  --num-patterns|n <Integer>
    # Total number of pairs pattern-text to generate.

  --length|l <Integer>
    # Length of the generated pattern (ie., query or sequence)
    #  and text (i.e., target or reference)

  --pattern-length|P <Integer>
    # Length of the generated pattern (ie., query or sequence)

  --text-length|T <Integer>
    # Length of the generated text (i.e., target or reference)

  --error|e <Float>
    # Total error-rate between the pattern and the text
    #  (allowing single-base mismatches, insertions and deletions).
    #  This parameter may modify the final length of the text.

  --debug|g
    # Output debug information.

  --help|h
    # Output a succinct manual for the tool.
```

## ALIGNMENT BENCHMARK TOOL

### Introduction to Alignment Benchmarking

The *align-benchmark* tool is provided to test and compare the performance of various pairwise alignment implementations.
This tool takes as input a dataset containing pairs of sequences (i.e., pattern and text) to align.
Patterns are preceded by the '>' symbol and texts by the '<' symbol. Example:

```text
>ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
<GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
>CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
<CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC
[...]
```

Once you have the dataset ready, you can run the *align-benchmark* tool to benchmark the performance of a specific pairwise alignment method.
For example, the QuickEd algorithm:

``` bash
$> ./bin/align_benchmark -i sample.dataset.seq -a quicked
...processed 100000 reads (alignment = 76431.865 seq/s)
...processed 200000 reads (alignment = 76604.137 seq/s)
[...]
...processed 5000000 reads (alignment = 69243.548 seq/s)
[Benchmark]
=> Total.reads              5000000
=> Time.Benchmark           1.20 m  (    1   call,  72.21  s/call {min72.21s,Max72.21s})
  => Time.Alignment        51.17 s  ( 70.86 %) (    5 Mcalls,  10.23 us/call {min8.80us,Max222.98us})
```

The *align-benchmark* tool will finish and report the overall benchmark time (including reading the input, setup, checking, etc.) and the time taken by the algorithm (i.e., *Time.Alignment*).
If you want to measure the accuracy of the alignment method, you can add the option `--check` with some of the available options and all the alignments will be verified.

```bash
$> ./bin/align_benchmark -i sample.dataset.seq -a quicked --check "alignment"
...processed 100000 reads (alignment = 21173.175 seq/s)
...processed 200000 reads (alignment = 20512.330 seq/s)
[...]
...processed 5000000 reads (alignment = 17404.361 seq/s)
[Benchmark]
=> Total.reads              5000000
=> Time.Benchmark           4.79 m  (    1   call, 287.28  s/call {min287.28s,Max287.28s})
  => Time.Alignment         1.01 m  ( 21.02 %) (    5 Mcalls,  12.08 us/call {min9.24us,Max630.13us})
[Accuracy]
 => Alignments.Correct        5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
 => Score.Correct             5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => Score.Total            24.20 Mscore uds.            (samples=5M{mean4.84,min0.00,Max5.00,Var0.00,StdDev0.00)}
     => Score.Diff            0.00 score uds.  (  0.00 %) (samples=0,--n/a--)}
 => CIGAR.Correct             4.75 Malg        ( 94.92 %) (samples=4M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
 => CIGAR.Breakdown
   => CIGAR.Matches         483.76 Mbases      ( 96.75 %) (samples=483M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Mismatches        8.29 Mbases      (  1.66 %) (samples=8M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Insertions        7.95 Mbases      (  1.59 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Deletions         7.96 Mbases      (  1.59 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
```

If you want even more detailed information about the alignment (e.g time breakdown of algorithm steps), you can add the option `--verbose`:

```bash
$> ./bin/align_benchmark -i sample.dataset.seq -a quicked --verbose
...processed 100000 reads (alignment = 75552.554 seq/s)
...processed 200000 reads (alignment = 75619.252 seq/s)
[...]
...processed 5000000 reads (alignment = 59000.338 seq/s)
[Benchmark]
=> Total.reads              5000000
=> Time.Benchmark           1.41 m  (    1   call,  84.75  s/call {min84.75s,Max84.75s})
  => Time.Alignment         1.01 m  ( 71.67 %) (    5 Mcalls,  12.15 us/call {min9.01us,Max961.50us})
  => Time.Windowed Small   26.38 s  ( 31.12 %) (    5 Mcalls,   5.28 us/call {min3.79us,Max942.88us})
  => Time.Windowed Large       0 ns (  0.00 %) (    0  calls,   n/a   s/call)
  => Time.Banded               0 ns (  0.00 %) (    0  calls,   n/a   s/call)
  => Time.Align            31.24 s  ( 36.87 %) (    5 Mcalls,   6.25 us/call {min4.64us,Max822.08us})
```

> [!NOTE]
> The overall benchmark time will increase due to the overhead introduced by the checking routine, however the *Time.Alignment* should remain the same.

### Algorithms & Implementations Summary

|          Algorithm Name           |       Code-name       |      Library      |
|-----------------------------------|:---------------------:|:-----------------:|
| **QuickEd**                       | quicked               | QuickEd           |
| **BandEd**                        | edit-banded           | QuickEd           |
| **WindowEd**                      | edit-windowed         | QuickEd           |
| **Edlib**                         | edlib                 | Edlib             |
| **DP Edit**                       | edit-dp               | None *(in-place)* |
| **DP Edit (Banded)**              | edit-dp-banded        | None *(in-place)* |
| **Bit-Parallel-Myers (BPM)**      | edit-bpm              | None *(in-place)* |

> *DP*: Dynamic Programming

### Command-line Options (align-benchmark)

#### Input

```bash
  --algorithm|a <algorithm-code-name>
    # Selects pair-wise alignment algorithm/implementation.

  --input|i <File>
    # Filename/path to the input SEQ file. That is, a file containing the sequence pairs to
    #  align. Sequences are stored one per line, grouped by pairs where the pattern is
    #  preceded by '>' and text by '<'.

  --output|o <File>
    # Filename/path of the output file containing a brief report of the alignment. Each line
    #  corresponds to the alignment of one input pair with the following tab-separated fields:
    #  <SCORE>  <CIGAR>

  --output-full <File>
    # Filename/path of the output file containing a report of the alignment. Each line
    #  corresponds to the alignment of one input pair with the following tab-separated fields:
    #  <PATTERN-LENGTH>  <TEXT-LENGTH>  <SCORE>  <PATTERN>  <TEXT>  <CIGAR>
```

#### Algorithm parameters

For a detailed list of parameters for each algorithm, please refer to the main [README](../README.md#alignment-methods-and-parameters-inside-the-quicked-library).

```bash
  --bandwidth <INT>
    # Selects the bandwidth to use in the BandEd algorithm.

  --window-size <INT>
    # Select the window size to use in the WindowEd algorithm.

  --overlap-size <INT>
    # Select the overlap size to use in the WindowEd algorithm.

  --hew-threshold <INT>
    # Error percentage threshold for HEWs in QuickEd.
    #  This parameter sets the same threshold for both steps.

  --hew-percentage <INT>
    # Percentage of HEW to consider that the estimation is not fitted in QuickEd.
    #  This parameter sets the same percentage for both steps.

  --only-score
    # Only compute the alignment score.

  --force-scalar
    # Force the use of the scalar implementations of the algorithms.
```

#### Others

```bash
  --num-threads|t <INT>
    # Sets the number of threads to use to align the sequences. If the multithreaded mode is
    #  enabled, the input is read in batches of multiple sequences and aligned using the number
    #  of threads configured.

  --batch-size <INT>
    # Select the number of pairs of sequences to read per batch in the multithreaded mode.

  --check|c 'correct'|'score'|'alignment'
    # Activates the verification of the alignment results.

  --check|c 'display'
    # Activates the display of the alignment results as part of the printed output.

  --progress|P <INT>
    # Sets the number of sequences to process before printing a progress message.

  --verbose|v <INT>
    # Sets the verbosity level of the tool.

  --help|h
    # Output a succinct manual for the tool.
```

