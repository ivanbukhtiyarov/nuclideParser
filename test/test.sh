MAIN_DIR = $pwd
echo 'my dir $MAIN_DIR'
git clone https://github.com/google/googletest.git ./lib/googletest
mkdir cmake-build-debug && cd cmake-build-debug
cmake ..
make all

