clear;clc;
problem.domain=[0 2];
problem.f=1;
s=sfwgd1(problem,0);
s=s.solve(2);
s.lift
x=linspace(0,1);
plot(x,polyval(s.lift(1,:),x),x+1,polyval(s.lift(1,:),x+1))



