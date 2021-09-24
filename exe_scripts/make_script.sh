cd ../build
rm -r bin/ && rm CMakeCache.txt && rm -r CMakeFiles/ && rm -r src/ && rm Makefile && rm cmake_install.cmake
cmake .. -DNETCDF_INCLUDE_DIR=/mnt/localscratch/include/ -DCODIPACK_INCLUDEDIR=/mnt/localscratch/include/ -DCMAKE_BUILD_TYPE=release -DTARGET=simulation && make -j4
# -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DTRACE_QV:BOOL=ON
# -DTRACE_QV:BOOL=ON -DTRACE_QR:BOOL=ON -DTRACE_SAT:BOOL=ON -DTRACE_TIME:BOOL=ON
# -DTARGET=timing