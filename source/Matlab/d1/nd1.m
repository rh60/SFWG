clear;clc;
t = linspace(0,2,3);
h = t(2)-t(1);
%B = [h/2 1; h/3 h/2] % difference matrix
B = powerMatrix(1,h);
D = B^(-1);

rhs = @(b,u) [ diff(b); b(2:end)-u];


B2 = rhs([0 1 0],[0 0]);
B1 = rhs([0 0 0],[1 0]);
B3 = rhs([0 0 0],[0 1]);
w = D*[B1 B2 B3]

K = zeros(3);
for i=1:3
    for j=1:3
        K(i,j) = prod(w(:,2*i-1:2*i),w(:,2*j-1:2*j));
    end
end
b = zeros(3,1);
b(2) = prod([1 1],[0 0]);
b(1) = prod([1 1],[1 0]);
b(3) = prod([1 1],[0 1]);


sol = K\b


function r=prod(p,q)
P1=polyint(conv(p(:,1)',q(:,1)'));
P2=polyint(conv(p(:,2)',q(:,2)'));
r=polyval(P1,1)-polyval(P1,0)+polyval(P2,1)-polyval(P2,0);
end

function P = powerMatrix(k,h)
P=zeros(k+1);
p=1;
for i=1:k+1
    q=1;
    for j=1:i
        pint=polyint(conv(p,q));
        P(i,j)=polyval(pint,h)-polyval(pint,0);
        q=[q 0];
    end
    p=[p 0];
end
P=P+tril(P,-1)';
P=flip(P,2);
end

function p = lifting(sol)
n = numel(sol);
p = 0;
end