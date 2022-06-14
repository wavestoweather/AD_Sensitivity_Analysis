cd ../build
rm -r bin/ && rm CMakeCache.txt && rm -r CMakeFiles/ && rm -r src/ && rm Makefile && rm cmake_install.cmake
cmake .. -DNETCDF_INCLUDE_DIR=/mnt/localscratch/include/ -DCODIPACK_INCLUDEDIR=/mnt/localscratch/include/ -DCMAKE_BUILD_TYPE=release -DTARGET=simulation -DHDF5_DIR=/data/Downloads/hdf5/CMake-hdf5-1.12.1/build/HDF5-1.12.1-Linux/HDF_Group/HDF5/1.12.1 -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON && make -j4

# -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON  -DB_EIGHT:BOOL=ON -DDEVELOP:BOOL=ON -DTRACE_COMM:BOOL=ON -DTRACE_READER:BOOL=ON