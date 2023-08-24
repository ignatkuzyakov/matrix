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
