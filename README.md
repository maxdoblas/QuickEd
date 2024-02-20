# QuickEd

Fast pairwise alignment using bit-parallel edit-distance algorithms

## Build

```bash
mkdir build && cd build
```

```bash
cmake ..
make
```

## Testing

To execute all tests, run:

```bash
cd build
make test # (or ctest)
```

## Debugging

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

### Valgrind

Make sure that Valgrind is installed on your system.
Then run:

```bash
cd build
ctest -T memcheck
```

**Note:** Valgrind runs will fail for any test that uses SSE4.1 intrinsics (i.e. any test that uses either WINDOWED 2x1 or QUICKED, without forceScalar).
This is expected behaviour.

### AddressSanitizer and UndefinedBehaviorSanitizer

To enable ASAN and UBSAN, build using:

```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DASAN=ON ..
make
```

**Note:** ASAN and UBSAN runs will fail for any test that uses SSE4.1 intrinsics (i.e. any test that uses either WINDOWED 2x1 or QUICKED, without forceScalar) or the Python binding.
This is expected behaviour.
