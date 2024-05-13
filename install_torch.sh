LIB_TORCH=libtorch
TORCH_SPARSE=pytorch_sparse
TORCH_SCATTER=pytorch_scatter

curl https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.3.0%2Bcpu.zip > torch.zip
unzip torch.zip 
rm -r torch.zip
echo "----- installing torch scatter -----"
git clone https://github.com/rusty1s/pytorch_scatter.git $TORCH_SCATTER --recurse-submodules
cd $TORCH_SCATTER && mkdir build --parents && cd build
sed -i "s/CMAKE_CXX_STANDARD [1-2][0-9]/CMAKE_CXX_STANDARD 17/" ../CMakeLists.txt
echo $(realpath ../../$LIB_TORCH) 
echo $(realpath $(pwd)/../install)
cmake .. -DCMAKE_PREFIX_PATH="$(realpath ../../$LIB_TORCH)" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH="$(realpath $(pwd)/../install)"
make -j $(nproc)
make install
cd ../..

echo "----- installing torch sparse -----"
git clone https://github.com/rusty1s/pytorch_sparse.git $TORCH_SPARSE --recurse-submodules
cd $TORCH_SPARSE && mkdir build --parents && cd build
sed -i "s/CMAKE_CXX_STANDARD [1-2][0-9]/CMAKE_CXX_STANDARD 17/" ../CMakeLists.txt
echo $(realpath ../../$LIB_TORCH) 
echo $(realpath $(pwd)/../install)
cmake .. -DCMAKE_PREFIX_PATH="$(realpath ../../$LIB_TORCH)" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH="$(realpath $(pwd)/../install)"
make -j $(nproc)
make install

