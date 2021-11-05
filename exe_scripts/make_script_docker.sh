if [ ! -d ../build ]
then
    mkdir ../build
fi
cd ../build
ls /usr/lib/x86_64-linux-gnu/
cmake .. -DCMAKE_BUILD_TYPE=release -DTARGET=simulation  -DTRUSTED_DATA:BOOL=ON -DB_EIGHT:BOOL=ON && make -j4
