clc;

disp('Compilation via MEX of AMI source code written in Cpp');
disp('Compiling..');
mex AMIcpp.cpp
disp('Done.');
