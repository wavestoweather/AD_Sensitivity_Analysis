if [ ! -d ../build ]
then
    mkdir ../build
fi
cd ../build
cmake .. -DCMAKE_BUILD_TYPE=release -DTARGET=simulation  -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON && make -j4
