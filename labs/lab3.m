function [u] = stepseq(n0,nstart,nend)
n = nstart:nend;
u = [1*((n-n0)>=0)];
end

function delta = impseq(n0,nstart,nend)
n = nstart:nend;
delta = [1*((n-n0)==0)];
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
% magX=abs(sin(N*w/2)./sin(w/2));
% plot(w/pi,magX); grid

x=ones(1,N);
X=fft(x,N);
hold on;

stem(2*[0:N-1]/N,abs(X),'r');
x2=[ones(1,N) zeros(1,N)];
X2=fft(x2,2*N);
stem(2*[0:2*N-1]/(2*N),abs(X2),'k'); % 2*N equispaced samples

x6=[ones(1,N) zeros(1,5*N)];
X6=fft(x6,6*N);
stem(2*[0:6*N-1]/(6*N),abs(X6),'m'); % 6*N equispaced samples
hold off