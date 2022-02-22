% mex glmnet
fprintf('Compiling MEX files...\n');
mex -largeArrayDims -v glmnetMex.f GLMnet.f