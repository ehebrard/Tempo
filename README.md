# Tempo

Requires cmake

On the cluster, I am using cmake -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_BUILD_TYPE=release ..
On my Mac: cmake -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX14.4.sdk -DCMAKE_BUILD_TYPE=release ..

Then 'make testsearch' to build the main exec (make will build a couple more executables)