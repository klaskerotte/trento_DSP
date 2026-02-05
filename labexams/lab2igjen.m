N = 400;
n = 0:N-1;

x1 = 1/50 *sin(pi/25 *n) .* cos(pi/50 *n) .* sin(pi/2 * n);
x2 = 2*cos(4/5 * pi * n).*x1;

wfft = linspace(0,2*pi,512);
X1 = fft(x1,512);
X2 = fft(x2,512);

M = 500;
wt1 = 0.4;
wt2 = 0.6;
d = 0.01;
bands = [0 wt1-d wt1+d wt2-d wt2+d 1];
amps = [1 1 0 0 1 1];
Rp = 3;
Rs = 65;

[bpm,apm] = firpm(M,bands,amps);
[BPM,w] = freqz(b,a,512);

y1 = filter(bpm,apm,x1);
y2 = filter(bpm,apm,x2);