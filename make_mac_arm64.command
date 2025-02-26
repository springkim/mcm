#!/bin/zsh

cd "$(dirname "$0")"

clang++ -DNDEBUG -O3 -fomit-frame-pointer -std=c++11 -D_FILE_OFFSET_BITS=64 -o mcm Archive.cpp Huffman.cpp MCM.cpp Memory.cpp Util.cpp Compressor.cpp File.cpp LZ.cpp Tests.cpp -lpthread

echo "Build complete!"
sleep 1