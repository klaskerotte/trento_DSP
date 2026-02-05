function [delta, n] = impseq(n0,nstart,nend)
n = nstart:nend;
delta = [(n==n0)];
end


N = 200;
n = 0:N-1;

x = zeros(1,N);
for r = 1:25
    x = x + impseq(8*r-1,0,N-1);
end
X = fft(x,512);

w_c = 0.62;

M = 100;
[b,a] = fir1(M,w_c,"high");
[B,w] = freqz(b,a,512);

y = filter(b,a,x);

d = 0.05;
[order, wn] = buttord(w_c-d,w_c+d,3,65);
[bf,ba] = butter(order,wn,"high");

[BF,w] = freqz(bf,ba,512);

ybf = filter(b,a,x);

