cd ../build
rm -r bin/ && rm CMakeCache.txt && rm -r CMakeFiles/ && rm -r src/ && rm Makefile && rm cmake_install.cmake
cmake .. -DTRUSTED_DATA:BOOL=ON -DTARGET=simulation && make -j 6
