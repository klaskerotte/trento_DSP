function [u] = stepseq(n0,nstart,nend)
n = nstart:nend;
u = [1*((n-n0>=0))];
end

function [x,nx] = sigfold(a,na)
x = fliplr(a);
nx = -fliplr(na);
end


function [y,m] = sigshift(x,n,k)
m = n+k;
y = x;
end

function [delta] = impseq(n0,nstart,nend)
n = nstart:nend;
delta = [1*((n-n0)==0)];
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

function y = myconv(x,h)
% y=myconv(x,h) direct computation of the convolution between x and h
L=length(x)+length(h)-1; % output length
y=zeros(1,L); % output sequence
x(L)=0; %zero padding of x till length L
h(L)=0; %zero padding of h till length L
for i=1:L
    for j=1:i
        y(i)=y(i)+x(j)*h(i-j+1); % output computation
        % Beware of the indexes (starting from 1 and not from 0 !)
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1

% Given the equation y[n] =x[n] −0.2y[n − 1] + 0.6y[n − 2]
% y[-1] = 2, y[-2] = -1
% Compute y[n] for 0 =< n =< 50
n = 0:50;

a = [1 0.2 -0.6];
b = [1];

x = 0.7.^n.*stepseq(0,0,50);

ic = filtic(b,a,[2,-1]);

y = filter(b,a,x,ic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2

%y[n]+1.12y[n-1]=0.1x[n]+0.2x[n-1], compute the zero input response (free response).

a = [1 1.12];
b = [0.1 0.2];

x = 0*stepseq(0,0,50);

ic = filtic(b,a,[1]);
y = filter(b,a,x,ic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1

x = [3, 11, 7, 0, -1, 4, 2];
nx = [-3:3];

h = [2, 3, 0, -5, 2, 1];
nh = [-1:4];

[y, ny] = conv_m(x,nx,h,nh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2

%y[n] − y[n − 1] + 0.9y[n − 2] = x[n] +x[n − 1];

% A

n = -20:100;

a = [1 -1 0.9];
b = [1 1];

h = impz(b,a,n);
stableh = sum(abs(h));


% B
x = stepseq(0,-20,100);
s = filter(b,a,x);
stables = sum(abs(s));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3

% 𝒚[𝒏] − 𝟎.𝟓𝒚[𝒏 − 𝟏] + 𝟎.𝟕𝟓𝒚[𝒏 − 𝟐] = 𝟐.𝟓𝒙[𝒏] + 𝟑𝒙[𝒏 − 𝟏] + 𝟐𝒙[𝒏 − 𝟐]

% x_1[n]=cos(2πn/10) and x2[n]=sin(4πn/5)

n = 0:40;

a = [1 -0.5 0.75];
b = [2.5 3 2];

x1 = 0.1*cos(2*pi*n/10);
x2 = cos(4*pi*n/5);

ypast1 = [0 0];
ypast2 = [-1 2];

ic1 = filtic(b,a,ypast1);
ic2 = filtic(b,a,ypast2);
y1 = filter(b,a,x1,ic1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 5

x = [3 11 7 0 -1 4 2];
nx = -3:3;

[y,ny] = sigshift(x,nx,2);
w = randn(size(y));

y = y+w;


[x,nx] = sigfold(x,nx);

[rxy,nrxy] = conv_m(x,nx,y,ny);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6

% y[n] = 5y[n − 1] + x[n], y[−1] = 0, for x[n] = u[n], 0 ≤ n ≤ 100

n = 0:100;

a = [1 -5];
b = [1];

ypast = 0;
x = stepseq(0,0,100);


ystep = filter(b,a,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 7

% y[n] = y[n − 1] + y[n − 2] + x[n], y[−1] = y[−2] = 0

% A
n = 0:100;
a = [1 -1 -1];
b = [1];

ypast = [0 0];
x = impseq(0,0,100);

impresponse = filter(b,a,x,ypast);
% B the system is unstable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 8

%y[n] = 1.15y[n − 1] − 1.5y[n − 2] + 0.7y[n − 3] − 0.25y[n − 4] + 0.18x[n] + 0.1x[n − 1] + 0.3x[n − 2] + 0.1x[n − 3] + 0.18x[n − 4]

% A

n = 0:100;
b = [0.18 0.1 0.3 0.1 0.18];
a = [1 -1.15 1.5 -0.7 0.25];

% Part (a):
h = impz(b,a,length(n));

% Part (b):
u = stepseq(0,0,n(end));
y = filter(b,a,u);

% Part (c):
nconv = 0:length(u)+length(h)-2;
y = conv(h,u);

% Part (d):
y = filter(h,1,u);
stem(n,y)