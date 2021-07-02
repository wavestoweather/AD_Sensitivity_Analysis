cd ../build
rm -r bin/ && rm CMakeCache.txt && rm -r CMakeFiles/ && rm -r src/ && rm Makefile && rm cmake_install.cmake && cmake .. -DNETCDF_INCLUDE_DIR=/mnt/localscratch/include/ -DCODIPACK_INCLUDEDIR=/mnt/localscratch/include/ -DNPROCS=4 -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON && make -j4
