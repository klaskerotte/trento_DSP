%%MatLab journey. Lab 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta]= impseq(n0,nstart,nend)
n = nstart:nend;
delta = [1*[(n-n0)==0]];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u] = stepseq(n0,nstart,nend)
n = nstart:nend;
u = [1*[(n-n0 >= 0)]];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,m] = sigshift(x,n,k)
m = n + k;
y = x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,n] = sigadd(x1,n1,x2,n2)

% Determine full index range
n_start = min(n1(1), n2(1));
n_end   = max(n1(end), n2(end));
n = n_start:n_end;

% Pre-allocate
y = zeros(1, length(n));

% ---- Map x1 into full range ----
idx1_start = n1(1) - n_start + 1;
idx1_end   = idx1_start + length(x1) - 1;
y(idx1_start:idx1_end) = y(idx1_start:idx1_end) + x1;

% ---- Map x2 into full range ----
idx2_start = n2(1) - n_start + 1;
idx2_end   = idx2_start + length(x2) - 1;
y(idx2_start:idx2_end) = y(idx2_start:idx2_end) + x2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,n] = sigmult(x1,n1,x2,n2)

% Determine full index range
n_start = min(n1(1), n2(1));
n_end   = max(n1(end), n2(end));
n = n_start:n_end;

% Pre-allocate
y = zeros(1, length(n));

% Zero-padded versions of x1 and x2
x1_pad = zeros(1, length(n));
x2_pad = zeros(1, length(n));

% Map x1
idx1_start = n1(1) - n_start + 1;
idx1_end   = idx1_start + length(x1) - 1;
x1_pad(idx1_start:idx1_end) = x1;

% Map x2
idx2_start = n2(1) - n_start + 1;
idx2_end   = idx2_start + length(x2) - 1;
x2_pad(idx2_start:idx2_end) = x2;

% Pointwise multiplication
y = x1_pad .* x2_pad;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,n] = sigfold(x,n)
% implements y(n) = x(-n)
% -----------------------
% [y,n] = sigfold(x,n)
%
y = fliplr(x); n = -fliplr(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xe, xo, m] = evenodd(x,n)
% Real signal decomposition into even and odd parts
% -------------------------------------------------
% [xe, xo, m] = evenodd(x,n)
%
if any(imag(x) ~= 0)
    error('x is not a real sequence')
end
m = -fliplr(n);
m1 = min([m,n]); m2 = max([m,n]); m = m1:m2;
nm = n(1)-m(1); n1 = 1:length(n);
x1 = zeros(1,length(m));
x1(n1+nm) = x; x = x1;
xe = 0.5*(x + fliplr(x)); xo = 0.5*(x - fliplr(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1
n_e1 = -5:5;
x_e1 = impseq(-2,-5,5) + 2*impseq(3,-5,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2
n_e2 = 0:20;
x_e2 = n_e2.*(stepseq(0,0,20)-stepseq(10,0,20)) + 10.*exp(-0.3.*(n_e2-10)).*(stepseq(10,0,20)-stepseq(20,0,20));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3
n_e3 = 0:200;
x_e3 = 3*cos(0.04*pi.*n_e3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 4

n_e4 = -30:29;
x_e4 = (0:9) + 1;

ones_base = ones(10,6);
joint = (ones_base'.*x_e4)';
unsqueezed = joint(:)';
xthilde_e4 = unsqueezed;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 5
x_e5 = [1:7 6:-1:1];
n_e5 = -2:(length(x_e5)-2);

[x_e5_a1,n_e5_a1] = sigshift(x_e5,n_e5,5);
[x_e5_a2,n_e5_a2] = sigshift(x_e5,n_e5,-4);
[A5,nA5] = sigadd(2.*x_e5_a1,n_e5_a1,-3.*x_e5_a2,n_e5_a2);
%stem(nA5,A5);

[x_e5_b1,n_e5_b1] = sigfold(x_e5,n_e5);
[x_e5_b1,n_e5_b1] = sigshift(x_e5_b1,n_e5_b1,3);

[x_e5_b2,n_e5_b2] = sigshift(x_e5,n_e5,2);
[x_e5_b2,n_e5_b2] = sigmult(x_e5_b2,n_e5_b2,x_e5,n_e5);

[B5,nB5] = sigadd(x_e5_b1, n_e5_b1, x_e5_b2,n_e5_b2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6

n_e6 = -10:10;
x_e6 = zeros(1,length(n_e6));

for k = -5:5
    x_e6 = x_e6+ exp(-abs(k)).*(impseq(2*k,-10,10));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 7

n_e7 = -30:30;
x_e7 = exp((-0.05+0.3j)*n);


% subplot(2,2,1);
% stem(n_e7,real(x_e7));
% subplot(2,2,2);
% stem(n_e7,imag(x_e7));
% subplot(2,2,3);
% stem(n_e7,abs(x_e7));
% subplot(2,2,4);
% stem(n_e7,angle(x_e7));

%polar(angle(x_e7),abs(x_e7),'-o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 8

n_e8 = -10:20;
x_e8 = stepseq(0,-10,20) - stepseq(10,-10,20);

[x_e8_even,x_e8_odd,n_e8]= evenodd(x_e8,n_e8);

%subplot(2,1,1)
%stem(n_e8,x_e8_even)

%subplot(2,1,2)
%stem(n_e8,x_e8_odd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 9

n1 = -5:15;
x1 = 3.*impseq(-2,-5,15) + 2.*impseq(0,-5,15) - impseq(3,-5,15) + 5.*impseq(7,-5,15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2 = -5:15;
x2 = 10.*stepseq(0,-5,15) - 5*stepseq(5,-5,15) - 10.*stepseq(10,-5,15) + 5.*stepseq(15,-5,15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n3 = -20:20;
x3_1 = exp(0.1.*n3);
x3_2 = (stepseq(-20,-20,20)-stepseq(10,-20,20));
x3 = sigmult(x3_1,n3,x3_2,n3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4 = -200:200;
x4 = 2*(cos(0.49*pi*n4) + cos(0.51*pi* n4));
% stem(n4,x4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n5 = n4;
x5 =  2*sin(0.01*pi*n5) .* cos(0.5*pi*n5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n6 = 0:100;
x6 = exp(-0.05*n6).*sin(0.1*pi.*n6 + pi/3);
% stem(n6,x6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x7 = [3 1 0 -1 2];
periods = 10;
basex7 = ones(length(x7),periods);
x7 = (x7'.*basex7);
x7 = x7(:)';
n7 = 0:length(x7)-1;
% stem(n7,x7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 10
x = [2 4 -3 1 -5 4 7]; % n = 0 at matlab-index 4 (1)
xn = [-3:3];

[x11,n11] = sigshift(x,xn,3);
[x12,n12] = sigshift(x,xn,-4);
[x1,n1] = sigadd(2*x11,n11,3*x12,n12);
[x1,n1] = sigadd(x1,n1,-x,xn);
% stem(n1,x1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x211,n211] = sigshift(x,xn,-3);
[x212,n212] = sigshift(x,xn,2);

[x221,n221] = sigfold(x,xn);
[x221,n221] = sigshift(x221,n221,1);
[x222,n222] = sigshift(x,xn,-1);

[x21,n21] = sigmult(x211,n211,x212,n212);
[x22,n22] = sigmult(x221,n221,x222,n222);

[x2,n2] = sigadd(x21,n21,x22,n22);
% stem(n2,x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n3 = -10:10;

x311 = 2*exp(0.5*n3);

x312 = x;
n312 = xn;

[x31,n31] = sigmult(x311,n3,x312,n312);

x321 = cos(0.1*pi*n3);
n321 = n3;
[x322,n322] = sigshift(x,xn,-2);

[x32,n32] = sigmult(x322,n322,x321,n321);

[x,xn] = sigadd(x31,n31,x32,n32);
% stem(xn,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 11

x1 = 0:9;
n1= 0:9;

[x1e,x1o,n1] = evenodd(x1,n1);

%subplot(2,2,1)
%stem(n1,x1e)
%subplot(2,2,2)
%stem(n1,x1o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n2 = -10:10;

x21 =  exp(0.1*n);

x221 = stepseq(-5,-10,10);
x222 = stepseq(10,-10,10);

x22 = sigadd(x221,n2,x222,n2);

x2 = sigmult(x21,n2,x22,n2);

[x2e,x2o,n2] = evenodd(x2,n2);

% subplot(2,2,1)
% stem(n2,x2e)
% subplot(2,2,2)
% stem(n2,x2o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n3 = -20:20;
x3 = cos(0.2*pi*n3 + pi/4);

[x3e,x3o,n3] = evenodd(x3,n3);


% subplot(2,2,1)
% stem(n3,x3e)
% subplot(2,2,2)
% stem(n3,x3o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n4 = 0:100;

x41 = exp(-0.05*n4);
x421 = sin(0.1*pi*n4 + pi/3);
x422 = stepseq(0,0,100);
x42 = sigmult(x421,n4,x422,n4);

x4 = sigmult(x42,n4,x41,n4);

[x4e,x4o,n4] = evenodd(x4,n4);

subplot(2,2,1)
stem(n4,x4e)
subplot(2,2,2)
stem(n4,x4o);