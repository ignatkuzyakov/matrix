## Testing matrix


for google tests

```
sudo apt-get install libgtest-dev googletest googletest-tools google-mock
```

build

```
cmake -G Ninja -B build
cmake --build build
env CTEST_OUTPUT_ON_FAILURE=1 cmake --build build --target test
```

valgrind output

```
valgrind ./build/tests

==2562==
==2562== HEAP SUMMARY:
==2562==     in use at exit: 0 bytes in 0 blocks
==2562==   total heap usage: 1,124 allocs, 1,124 frees, 149,827 bytes allocated
==2562==
==2562== All heap blocks were freed -- no leaks are possible
==2562==
==2562== For lists of detected and suppressed errors, rerun with: -s
==2562== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```

