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
ctest -E edlib_tests --output-on-failure # Exclude Edlib tests
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

### AddressSanitizer and UndefinedBehaviorSanitizer

To enable ASAN and UBSAN, build using:

```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DASAN=ON ..
make
```

**Note:** ASAN and UBSAN could fail when executing the Python binding. This is known behaviour.