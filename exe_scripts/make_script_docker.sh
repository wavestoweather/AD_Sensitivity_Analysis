if [ ! -d ../build ]
then
    mkdir ../build
fi
cd ../build
cmake .. -DCMAKE_BUILD_TYPE=release -DTARGET=simulation  -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON && make -j4

if [ ! -d ../build_interface ]
then
    mkdir ../build_interface
fi
cd ../build_interface
cmake .. -DCMAKE_BUILD_TYPE=release -DTARGET=python_interface  -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON && make -j4
