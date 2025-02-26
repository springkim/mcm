mkdir WindowsX64_VS2019
cl /EHsc /O2 /DNDEBUG /DWIN32 /D_FILE_OFFSET_BITS=64 /MT /arch:SSE2 /Fe:mcm.exe Archive.cpp Huffman.cpp MCM.cpp Memory.cpp Util.cpp Compressor.cpp File.cpp LZ.cpp Tests.cpp
move mcm.exe .\WindowsX64_VS2019
DEL *.obj
pause