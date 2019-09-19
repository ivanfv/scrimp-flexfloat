cd flexfloat
rm -r build 2> /dev/null 
mkdir build
cd build
cmake -DENABLE_STATS=ON .. 
make
