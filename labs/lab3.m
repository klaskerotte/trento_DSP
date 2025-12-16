function [u] = stepseq(n0,nstart,nend)
n = nstart:nend;
u = [1*((n-n0)>=0)];
end

function delta = impseq(n0,nstart,nend)
n = nstart:nend;
delta = [1*((n-n0)==0)];
end

function [y] = fastconv(x,h)
% Fast convolution using FFT
% --------------------------------------------------------------
% [y] = fastconv(x,h)
% y = output sequence
% x = first input sequence
% h = second input sequence
%
L=length(x)+length(h)-1; % output length
N = 2^(ceil(log10(L)/log10(2))); % next power of 2
y = ifft(fft(h,N).*fft(x,N)); % zero padding up to length N
y=y(1:L);
end

function y = cirshift(x,m,N)
% Circular shift of m samples wrt size N
% in sequence x
% -----------------------------------------------------------
% [y] = cirshftt(x,m,N)
% y = output sequence containing the circular shift
% x = input sequence of length <= N
% m = shift in samples
% N = size of circular buffer
if length(x) > N
    error('N must be >= length of x')
end
x = [x zeros(1,N-length(x))];
nx = [0:N-1];
ny = mod(nx-m,N); % apply shift to sample indexes
y = x(ny+1); % indexes of Matlab vectors must start from 1
end

function y = circonvt(x1,x2,N)
% N-point circular convolution between x1 and x2: (time-domain)
% -------------------------------------------------------------
% [y] = circonvt(x1,x2,N)
% y = output sequence containing the circular convolution
% x1 = input sequence of length N1 <= N
% x2 = input sequence of length N2 <= N
% N = size of circular buffer
% Method: y(n) = sum (x1(m)*x2((n-m) mod N))
if length(x1) > N
    error('N must be >= the length of x1')
end
if length(x2) > N
    error('N must be >= the length of x2')
end
x1=[x1 zeros(1,N-length(x1))];
x2=[x2 zeros(1,N-length(x2))];
m = [0:1:N-1];
x2n = x2(mod(-m,N)+1); % obtain the sequence x2n[n]=x2[-n]
for n = 0:N-1
    y(n+1)=x1 * cirshftt(x2n,n,N).';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Lab 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1

%Determine the frequency response 𝑯 𝒆𝒋𝝎 of a system characterized by the impulse response
%h[n] = (0.9)^n u[n]
%Plot the magnitude and the phase responses.

% H[e^(j*omega)] = sum_0^inf(0.9*n*e^(j*omega*n)) --> Geometric series sum where |r| < 1  --> Analytical solution yields: H[e^(j*omega)] = 1/(1-0.9e^(j*omega))

omega = [0:500]*pi/500;
H = 1./(1-0.9*exp(1j*-omega));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 2

% X(e^{j*omega}) = sum_0^{N-1} e^{-j*omega*n} ----> Finite geometric series sum ----> X(e^(j*omega)) = (1-e^(-j*omega*N)) / (1-e^(-j*oemga))
% This gives |X(e^(j*omega)| =  abs(sin(omega*N/2) / sin(omega / 2))

N=11;
w=0:1/1000:2*pi;
magX=abs(sin(N*w/2)./sin(w/2));
% plot(w/pi,magX); grid

x=ones(1,N);
X=fft(x,N);
% hold on;
% stem(2*[0:N-1]/N,abs(X),'r');

x2=[ones(1,N) zeros(1,N)];
X2=fft(x2,2*N);
% stem(2*[0:2*N-1]/(2*N),abs(X2),'k'); % 2*N equispaced samples

x6=[ones(1,N) zeros(1,5*N)];
X6=fft(x6,6*N);
% stem(2*[0:6*N-1]/(6*N),abs(X6),'m'); % 6*N equispaced samples
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 3

n = -10:10;
x = (-0.9).^n;

Xf = fft(x,512);
omega = linspace(0,pi,512);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 4

n = -1:4;
x = [1 2 3 4 5];


x = circshift(x,-1);
n = 0:length(x)-1;

Xf = fft(x,512);
omega = linspace(0,2*pi,512);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 5

x1 = [1 1 1 1 1];
n1 = 0:4;

x2 = [1 1 1 1 1];
n2 = -2:2;

x2 = circshift(X2,-2);
n2 = 0:4;

X1f = fft(x1,512);
X2f = fft(x2,512);

omega = linspace(0,2*pi,512);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6

n = 0:32;
omega0 = 2*pi*5/32;
x = cos(omega0*n);

X1f32 = fft(x,32);
X1f64 = fft(x,64);
X1f128 = fft(x,128);

omega32 = linspace(0,2*pi,32);
omega64 = linspace(0,2*pi,64);
omega128 = linspace(0,2*pi,128);

% subplot(3,2,1)
% stem(omega32,abs(X1f32))
% subplot(3,2,3)
% stem(omega64,abs(X1f64))
% subplot(3,2,5)
% stem(omega128,abs(X1f128))

n = 0:32;
omega0 = 2*pi*5.5/32;
x = cos(omega0*n);

X1f32 = fft(x,32);
X1f64 = fft(x,64);
X1f128 = fft(x,128);

omega32 = linspace(0,2*pi,32);
omega64 = linspace(0,2*pi,64);
omega128 = linspace(0,2*pi,128);

% subplot(3,2,2)
% stem(omega32,abs(X1f32))
% subplot(3,2,4)
% stem(omega64,abs(X1f64))
% subplot(3,2,6)
% stem(omega128,abs(X1f128))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 7

[x,fs]=audioread('bluewhale.wav'); %read audio data from file

b= x(2.45e4:3.10e4);
B=fft(b);
f=fs/length(b)*[0:(length(b)-1)];
% plot(f/10,abs(B)) % plot with the true frequency scale

winlen=512; step=256; Nfft=1024;
% spectrogram(x,winlen,step,Nfft,fs/10,'yaxis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 8

n = -50:50;
x1 = ((0.7).^abs(n)).*(stepseq(-20,-50,50)-stepseq(21,-50,50));
x1 = circshift(x1,-20,512);
x1f = fft(x1,512);

x2 = n.*((0.9).^n).*(stepseq(0,-50,50)-stepseq(21,-50,50));
x2 = circshift(x2,-20,512);
x2f = fft(x2,512);

omega = linspace(0,2*pi,512);

plot(omega,abs(x1f))