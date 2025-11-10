function [u,n] = stepseq(n0,n1,n2)
n = n1:n2;
u = 1*[n-n0>=0];
end

function [d,n] = impseq(n0,n1,n2)
n = n1:n2;
d = 1*[n-n0==0];
end

function [y,ny] = conv_m(x,nx,h,nh)
% Modified convolution routine for signal processing
% --------------------------------------------------
% [y,ny] = conv_m(x,nx,h,nh)
% [y,ny] = convolution result
% [x,nx] = first signal
% [h,nh] = second signal
%
nyb = nx(1)+nh(1); nye = nx(length(x)) + nh(length(h));
ny = [nyb:nye];
y = conv(x,h);
end

% Difference equation is --> y[n] = x[n] -0.2y[n-1] + 0.6y[n-2]
% past history of y is --> y[-1] = 2, y[n-2] = -1
% input is x[n] = 0.7^n*u[n]

n = 0:50;

%define y-coefficients a (NB! remember to separate y and x on each side)
a = [1 0.2 -0.6];

%define x-coefficients b
b = [1];

%initialize past of y, i.e. y[-1] and y[-2]
ypast = [2 -1];

%Get initial filterstate
ic = filtic(b, a,ypast);

%Define inputsequence x[n]
x = 0.7.^n.*stepseq(0,0,50);

%solve diff-eq
y = filter(b,a,x,ic);



%New diff-eq y[n] + 1.12y[n-1] = 0.1x[n] + 0.2x[n-1]
a = [1 1.12];
b = [0.1 0.2];

ypast = [1]
ic = filtic(b,a,ypast);

y = filter(b,a,zeros(1,20),ic);



%Compute the convolution of x and h
x = [3, 11, 7, 0, -1, 4, 2];
nx = [-3:3];

h = [2, 3, 0, -5, 2, 1];
nh = [-1:4];

[y,ny] = conv_m(x,nx,h,nh);



%Given the following diff-eq --> y[n] − y[n − 1] + 0.9y[n − 2] = x[n] +x[n − 1]; ∀n
%A.
n = [-20:100];
a = [1, -1, 0.9];
b = [1 1];
h = impz(b,a,n);
subplot(2,1,1); stem(n,h);
title('Impulse Response'); xlabel('n'); ylabel('h(n)')

%B
n = [-20:100];
a = [1, -1, 0.9];
b = [1 1];
h = stepz(b,a,n);
subplot(2,1,1); stem(n,h);
title('Impulse Response'); xlabel('n'); ylabel('h(n)')

%C
