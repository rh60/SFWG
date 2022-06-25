clear;clc;
problem.domain=[0 2];
problem.f=1;
s=sfwgd1(problem,0);
s=s.solve(2);
s.sol

