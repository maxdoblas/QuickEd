QuickEd
&nbsp;
[![Release](https://img.shields.io/github/release/maxdoblas/quicked.svg)](https://github.com/maxdoblas/quicked/releases/latest)
[![Publication](https://img.shields.io/badge/Published%20in-BioRxiv-167DA4.svg)](https://doi.org/10.1101/2024.09.13.612714)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14526794.svg)](https://doi.org/10.5281/zenodo.14526794)
=====

QuickEd is a high-performance exact sequence alignment based on the bound-and-align paradigm.

Calculating the edit distance of two DNA sequences using the quicked library is as simple as:

```c
quicked_params_t params = quicked_default_params();
quicked_new(&aligner, &params);
quicked_align(&aligner, pattern, strlen(pattern), text, strlen(text));
```

QuickEd has bindings available for **C++** at [bindings/cpp](bindings/cpp) and **Python** at [bindings/python](bindings/python).

## Features

* Calculates **edit distance (Levenshtein distance)**.
* It can find **optimal alignment path** (instructions on how to transform the first sequence into the second sequence).
* Supports **multiple alignment algorithms**:
  * QuickEd
  * BandEd
  * WindowEd
* Efficient implementation of the bound-and-align paradigm.
* It scales to **long and noisy sequences**, even when finding an alignment path, while consuming very little memory thanks to the Hirschberg algorithm.

Contents
--------

* [Features](#features)
* [Contents](#contents)
* [Integrate QuickEd in your tools](#integrate-quicked-in-your-tools)
  * [Install QuickEd locally](#install-quicked-locally)
  * [Use QuickEd in your CMake project](#use-quicked-in-your-cmake-project)
  * [Manually add QuickEd files to your project](#manually-add-quicked-files-to-your-project)
* [Building](#building)
* [Usage](#usage)
  * [Configuring QuickEd](#configuring-quicked)
  * [Handling result of quicked\_align()](#handling-result-of-quicked_align)
  * [Alignment methods and parameters inside the QuickEd library](#alignment-methods-and-parameters-inside-the-quicked-library)
* [Testing](#testing)
* [Development and Debugging](#development-and-debugging)
  * [Generate code coverage data](#generate-code-coverage-data)
  * [AddressSanitizer and UndefinedBehaviorSanitizer](#addresssanitizer-and-undefinedbehaviorsanitizer)
* [Experimental Evaluation](#experimental-evaluation)

## Integrate QuickEd in your tools

You can use QuickEd in your project by linking the QuickEd static library into your binaries (check [Building](#building) for instructions on how to build the QuickEd static library).
Besides that, the only thing you have to do is to include `quicked.h` in your source files.

For examples on how to use QuickEd in your own code, look at the [examples](examples) folder, where we have several use case example codes for C, and for C++ and Python through the provided bindings.

The simplest implementation that uses QuickEd looks like this:

```c
#include <quicked.h>
#include <stdio.h>
#include <string.h>

int main(void) {
    quicked_aligner_t aligner;                          // Aligner object
    quicked_params_t params = quicked_default_params(); // Get a set of sensible default parameters.
    quicked_new(&aligner, &params);                     // Initialize the aligner with the given parameters

    const char* pattern = "ACGT";                       // Pattern sequence
    const char* text = "ACTT";                          // Text sequence
    quicked_align(&aligner, pattern, strlen(pattern), text, strlen(text));

    printf("Alignemnt: %d - $s\n", aligner.score, aligner.cigar);   // Print the score and CIGAR

    quicked_free(&aligner);                 // Free whatever memory the aligner allocated

    return 0;
}
```
### Install QuickEd locally

The easiest way to use QuickEd is to install the library files system-wide.

You can achieve this by running `sudo make install` after building.
Now, all you have to do to compile your project is `gcc aligner.c -o aligner -lquicked`.

> [!NOTE]
> You may need to have the necessary permissions to install the library on your system.

### Use QuickEd in your CMake project

Since QuickEd is a CMake project, you can easily include it in other CMake projects as a git submodule.

To do so, use `git submodule add https://github.com/maxdoblas/quicked external/quicked`
Then, you should add the following statements to your CMakeLists.txt file:

```cmake
add_subdirectory(external/quicked EXCLUDE_FROM_ALL)
target_link_libraries(your_exe quicked)
```

### Manually add QuickEd files to your project

> [!IMPORTANT]
> We advise against this method, as it is more error-prone and harder to maintain.
> You are on your own.

If the other methods don't suffice you, you can manually add the required QuickEd files to your proyect.

You need to copy the QuickEd static library `libquicked.a` and some header files to your project:

```text
your-project/
├─ quicked/
│  └─ quicked.h     <- Found at <this repo>/quicked/quicked.h
├─ quicked_utils/   <- Copy the indicated files from <this repo>/quicked_utils/include
│  └─ include/
│     ├─ commons.h
│     ├─ mm_allocator.h
│     ├─ profiler_timer.h
│     ├─ profiler_counter.h
│     └─ vector.h
├─ libquicked.a     <- Found at <this repo>/lib/libquicked.a
└─ aligner.c        <- Your code
```

Now you can compile it with `gcc aligner.c -o aligner -I quicked -I quicked_utils -L your-project -lquicked`.

## Building

To build QuickEd, clone the repository and update the submodules:

```bash
git clone https://github.com/maxdoblas/QuickEd.git -b benchmark
cd QuickEd
```

> [!IMPORTANT]
> QuickEd is built using CMake, so make sure you have it installed on your system. The minimum required version is 3.20

> [!NOTE]
> `benchmark` branch integrates [**A\* Pairwise Alignment**](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner). For more information about its building process, please refer to the original repo.

Execute the following command to build QuickEd using CMake:

```bash
cd tools/align_benchmark/external/astar-pairwise-aligner/
cargo build --release
cbindgen --lang c --cpp-compat --crate astarpa-c -o astarpa.h
cd ../../../../
mkdir build && cd build
cmake ..
make
```

This will create all the tools and example binaries in the `bin/` directory and the QuickEd static library (as well as the library objects for bindings) in the `lib/` directory.

Optionally, you can run.

```bash
sudo make install
```

To install the QuickEd library on your machine. This will install the QuickEd library and header files in your system's default location.

## Usage

The main function of QuickEd is `quicked_align`. Given two DNA strings (and their lengths), it will find the edit distance and alignment path.

```c
const char* pattern = "ACGT";                       // Pattern sequence
const char* text = "ACTT";                          // Text sequence
quicked_align(&aligner, pattern, strlen(pattern), text, strlen(text));
```

### Configuring QuickEd

`quicked_align` takes the `quicked_aligner_t` object. The `quicked_aligner_t` is created with the `quicked_new` constructor, which takes the configuration struct (`quicked_params_t`), allowing the customization of the different methods implemented inside the library. You can find all the information at [alignment methods](#alignment-methods-and-parameters-inside-the-quicked-library).

```c
typedef struct quicked_params_t {
    quicked_algo_t algo;
    unsigned int bandwidth;
    unsigned int window_size;
    unsigned int overlap_size;
    unsigned int hew_threshold[QUICKED_WINDOW_STAGES];
    unsigned int hew_percentage[QUICKED_WINDOW_STAGES];
    bool only_score;
    bool force_scalar;
    bool external_timer;
    mm_allocator_t *external_allocator;
} quicked_params_t;
```

You can easily use the `quicked_default_params()` function to create the default configuration. Then, you could set its different parameters members accordingly. For example, if we want to use the BandEd method with a 10% bandwidth, you have to do:

```c
quicked_aligner_t aligner;
quicked_params_t params = quicked_default_params();
params.algo = BANDED;  // Select the algorithm: Banded
params.bandwidth = 10; // Banded needs a bandwidth
quicked_new(&aligner, &params);
```

### Handling result of quicked_align()

The `quicked_align` function returns the status struct (`quicked_status_t`), and the alignment result (CIGAR and score) will be set inside the `quicked_aligner_t` object. We can check the status structure using the following functions:

```c
  quicked_status_t status = quicked_align(&aligner, pattern, pattern_length, text, text_length);

  if (quicked_check_error(status)){
    fprintf(stderr, "%s", quicked_status_msg(status));
    exit(1);
  } else {
    printf("Alignemnt: %d - $s\n", aligner.score, aligner.cigar);   // Print the score and CIGAR
  }
```

> [!CAUTION]
> It is important to free the `quicked_aligner_t` object. It dynamically allocates memory for the different fields inside when created. You can free it using the `quicked_free` function.

### Alignment methods and parameters inside the QuickEd library

The `quicked_params_t` configuration struct has the following parameters:

* **quicked_algo_t** `algo`: contains the method the aligning function uses. The methods implemented are:
  * `QUICKED`: The default method in the library. QuickEd is the edit distance implementation of the bound-and-align paradigm.
  * `WINDOWED`: The WindowEd algorithm is a fast and efficient heuristic algorithm that uses overlapping windows to compute the alignment.
  * `BANDED`: BandEd algorithm computes all the different cells with a score lower than a cutoff value provided with the `bandwidth` parameter.
  * `HIRSCHBERG`: It uses the BandEd algorithm with Hirschebrg optimization to reduce the memory footprint.

* **unsigned int** `bandwidth`: sets the bandwidth (or cutoff score) in % used on the BandEd algorithms.
* **unsigned int** `window_size`: sets the window size in blocks for the WindowEd algorithms. Also, it sets the window size for the WindowEd(L) inside the QuickEd method. **Note**: the size in cells will be `64*window_size`.
* **unsigned int** `overlap_size`: sets the overlap size in blocks for the WindowEd algorithms. Also, it sets the window size for the WindowEd(L) inside the QuickEd method. **Note**: the size in cells will be `64*window_size`.
* **unsigned int** `hew_threshold[2]`: The error percentage threshold inside a window to be considered a high error window (HEW). This parameter is only used inside Quicked. Position [0] refers to the WindowEd(S) step and position [1] to the WindowEd(L) step.
* **unsigned int** `hew_percentage[2]`: percentage of HEW in a particular WindowEd alignment to consider that the estimation is not fitted. This parameter is only used inside Quicked. Position [0] refers to the WindowEd(S) step and position [1] to the WindowEd(L) step.
* **bool** `force_scalar`: If set to true, it forces WindowEd and BandEd implementation to use the scalar code.

> [!WARNING]
> **Experimental Parameters**
>
> * **bool** `external_timer`: If set to true, it uses external timers and avoids the generation of a new allocator.
> * **mm_allocator_t** `*external_allocator`: If it is set with an external allocator, this allocator will be used instead of creating a new allocator object.

## Testing

You can execute all the tests found in the CI by running the following commands:

```bash
cd build
ctest -E edlib_tests --output-on-failure # Exclude Edlib tests
```

## Development and Debugging

> [!WARNING]
> Use these commands at your own risk :D

For debugging, build using:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

### Generate code coverage data

```bash
cd build
ctest -T coverage
```

### AddressSanitizer and UndefinedBehaviorSanitizer

To enable ASAN and UBSAN, build using:

```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DASAN=ON ..
make
```

> [!WARNING]
> There are some known issues with ASAN and UBSAN in general. Check https://stackoverflow.com/a/77895910 for more information.

> [!WARNING]
> ASAN and UBSAN could fail when executing the Python binding. This is known behavior.

## Experimental Evaluation

We evaluated the performance of QuickEd compared to other state-of-the-art sequence alignment libraries. In this section, we provide a brief overview of the steps to reproduce the experimental evaluation found in the paper.

### Datasets

The used datasets are enumerated in S1.2 of the Supplementary Material, and simulated datasets can be found at [Zenodo](https://doi.org/10.5281/zenodo.14526793).

### Methods

We evaluated QuickEd against the following libraries/algorithms:

* [Scrooge](https://github.com/CMU-SAFARI/Scrooge) [in _WFA2-lib_ [455be3b](https://github.com/smarco/WFA2-lib/tree/455be3bf8ba43f4fd7a33017b306a1bb5d8dcbec)]
* [Edlib](https://github.com/Martinsos/edlib) [[931be2b](https://github.com/Martinsos/edlib/tree/931be2b0909985551eb17d767694a6e64e31ebfa)]
* [KSW2](https://github.com/lh3/ksw2) [in _WFA2-lib_ [455be3b](https://github.com/smarco/WFA2-lib/tree/455be3bf8ba43f4fd7a33017b306a1bb5d8dcbec)]
* [A* Pairwise Aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) (_A*PA_ and _A*PA2_) [[8ca0f21](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner/tree/8ca0f216886765a85a2260e6ac0a6b2f3f69819e)]
* [WFA2-lib](https://github.com/smarco/WFA2-lib) (_WFA_ and _BiWFA_) [[42f8bf3](https://github.com/smarco/WFA2-lib/commit/42f8bf3178fb5fd5896f46ffdd08f6d5fb927b51)]

Edlib is integrated in the [**align_benchmark**](https://github.com/maxdoblas/QuickEd/tree/benchmark/tools) tool, as a ground-truth comparison.
A* Pairwise Aligners is also integrated through its C binding.

WFA/BiWFA can be used either in QuickEd's integration or through its original repository.

Scrooge and KSW2 can be found already integrated in the WFA2-lib repository (`benchmark` branch), and they have been evaluated using its [**align_benchmark**](https://github.com/smarco/WFA2-lib/tree/benchmark/tools) tool which is very similar to ours.

### Environment

All the experiments were run on a machine with the following specifications:

* **CPU**: Intel Xeon Platinum 8480+ @ 3GHZ _(fixed at this frequency)_
* **RAM**: 256 GB

All runs were performed using a single thread and a single, exclusive node (enforced by SLURM with `#SBATCH --hint=nomultithread` and `#SBATCH --exclusive`).

We also used the following versions for compilers and other tools:

* **GCC**: 12.3.0
* **CMake**: 3.25.1
* **Rust**: 1.85.0-nightly

With this setup, `AVX512` instructions are available whenever the tool can leverage them.
All libraries were compiled following the instructions provided in their respective repositories.

### Reproducing the results

To reproduce the results, you need to follow these steps:

**1.** Clone and compile both QuickEd and [WFA2-lib](https://github.com/smarco/WFA2-lib/tree/benchmark) repositories, according to the instructions provided in their README files, **using the `benchmark` branch**.

```bash
git clone -b 'benchmark' https://github.com/maxdoblas/QuickEd.git
cd QuickEd
cd tools/align_benchmark/external/astar-pairwise-aligner/
cargo build --release
cbindgen --lang c --cpp-compat --crate astarpa-c -o astarpa.h
cd ../../../../
mkdir build && cd build
cmake ..
make
```

```bash
git clone -b 'bechmark' https://github.com/smarco/WFA2-lib
cd WFA2-lib
make external-all
make all
```

**2.** Using the corresponding **align_benchmark** for the desired algorithm, run the following command:

```bash
./align_benchmark -i <dataset> -o <output> -a <algorithm> --verbose 4
```

The available algorithms in QuickEd repo are:

|       **Method**       | **`-a <algorithm>`** |
|:----------------------:|:--------------------:|
| QuickEd                | quicked              |
| * BandEd _(QuickEd)_   | edit-banded          |
| * WindowEd _(QuickEd)_ | edit-windowed        |
| WFA                    | wfa2-edit            |
| BiWFA                  | biwfa-edit           |
| A*PA                   | astarpa              |
| Edlib                  | edlib                |
| A*PA2 Simple           | astarpa2-simple      |
| A*PA2 Full             | astarpa2-full        |

Likewise, the available algorithms in WFA2-lib repo are:

|       **Method**      | **`-a <algorithm>`** |
|:---------------------:|:--------------------:|
| Scrooge               | scrooge              |
| KSW2 _(Gap-Affine)_   | ksw2-extz2-sse       |

For KSW2 to represent edit-distance, you need to add the `--affine-penalties 0,1,0,1` flag.

**3.** For _BandEd_ and _WindowEd_ methods, additional parameters can be set:

```bash
./align_benchmark -i <dataset> -o <output> -a edit-banded --bandwidth <int> --only-score
```

```bash
./align_benchmark -i <dataset> -o <output> -a edit-windowed -a edit-windowed --overlap-size <int> --window-size <int>
```

For _BandEd_ and _WindowEd_ executions without vectorization, you can use the `--force-scalar` flag.

**4.** To obtain the memory consumption, prepend the command with `/bin/time -v`:

```bash
/bin/time -v ./align_benchmark ...
```

Or alternatively with `\time`:

```bash
\time -v ./align_benchmark ...
```
