QuickEd
&nbsp;
[![Release](https://img.shields.io/github/release/mdoblas/quicked.svg)](https://github.com/mdoblas/quicked/releases/latest)
[![CI](https://img.shields.io/github/actions/workflow/status/mdoblas/quicked/build_and_test.yaml?branch=dev)](https://github.com/maxdoblas/QuickEd/actions/workflows/build_and_test.yaml)
[![Publication](https://img.shields.io/badge/Published%20in-BioRxiv-167DA4.svg)](https://www.biorxiv.org/)
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

* [](#)
  * [Features](#features)
  * [Contents](#contents)
  * [Using QuickEd in your project](#using-quicked-in-your-project)
    * [Approach #1: Using QuickEd header file and static library](#approach-1-using-quicked-header-file-and-static-library)
    * [Approach #2: Install the QuickEd library on the machine](#approach-2-install-the-quicked-library-on-the-machine)
    * [Approach #3: Use Quicked in your project via CMake as a git submodule](#approach-3-use-quicked-in-your-project-via-cmake-as-a-git-submodule)
  * [Building](#building)
  * [Usage and examples](#usage-and-examples)
    * [Configuring QuickEd](#configuring-quicked)
    * [Handling result of quicked\_align()](#handling-result-of-quicked_align)
  * [Alignment methods and parameters inside the QuickEd library](#alignment-methods-and-parameters-inside-the-quicked-library)
  * [Testing](#testing)
  * [Development and Debugging](#development-and-debugging)
    * [Generate code coverage data](#generate-code-coverage-data)
    * [AddressSanitizer and UndefinedBehaviorSanitizer](#addresssanitizer-and-undefinedbehaviorsanitizer)

## Using QuickEd in your project

You can use QuickEd in your project by linking the QuickEd library into your binaries (see [Building](#building) for instructions on how to build the QuickEd library).
In any case, the only thing that you have to do in your source files is to include `quicked.h`.

To use QuickEd in your project, look at our examples folder, where we have several use case example codes for C, C++, and Python.

Our simplest example has just one source file, the `basic.c` file, and it looks like this:

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

### Approach #1: Using QuickEd header file and static library

We need to copy the QuickEd static library and header files (check [Building](#building) on how to create a static library). We get the following project structure:

```text
quicked/         <- copied from quicked
├─ include/
│   └─ quicked.h
├─ quicked.a
└─ aligner.c     <- your program
```

Now you can compile it with `gcc aligner.c -o aligner -I quicked/include -L quicked -lquicked`.

### Approach #2: Install the QuickEd library on the machine

Alternatively, you could avoid copying any QuickEd files and instead install libraries by running `sudo make install` (check [Building](#building) for exact instructions). Now, all you have to do to compile your project is `gcc aligner.cpp -o aligner -lquicked`.

> [!NOTE]
> If you get an error message like `cannot open shared object file: No such file or directory`, ensure your linker includes the path where Quicked is installed.

### Approach #3: Use Quicked in your project via CMake as a git submodule

If you use CMake for infrastructure, you can include QuickEd as a git submodule with the command `git submodule adds https://github.com/mdoblas/quicked external/quicked`. Then, you should add the following statements to your CMakeLists.txt file:

```cmake
add_subdirectory(external/quicked EXCLUDE_FROM_ALL)
target_link_libraries(your_exe quicked)
```

## Building

QuickEd is built using CMake.

Execute the following command to build QuickEd using CMake:

```bash
mkdir build && cd build
cmake ..
make
```

This will create all the tools and example binaries in the `bin/` directory and the QuickEd libraries (static and shared) in the `lib/` directory.

Optionally, you can run.

```bash
sudo make install
```

To install the QuickEd library on your machine.

## Usage and examples

The main function in QuickEd is `quicked_align.` Given two sequences (and their lengths), it will find the edit distance and alignment path.

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

## Alignment methods and parameters inside the QuickEd library

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
* **bool** `only_score`: If set to true, turn off the CIGAR generation for the WindowEd and BandEd methods.
* **bool** `force_scalar`: If set to true, it forces the WindowEd implementation to use the scalar code.

> [!WARNING] Experimental Parameters
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
> ASAN and UBSAN could fail when executing the Python binding. This is known behavior.
