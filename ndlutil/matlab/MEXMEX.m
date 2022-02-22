% mex
fprintf('Compiling MEX files...\n');
mex -compatibleArrayDims lnDiffErfs.c
mex -compatibleArrayDims xgamrnd.c
%%mex -largeArrayDims minFunc/mcholC.c
%%mex -largeArrayDims minFunc/lbfgsC.c