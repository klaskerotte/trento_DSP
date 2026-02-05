%% Exercise 1

wc = 1;
N = 4;

[b,a] = butter(N,wc,'s');

[Bs,ws] = freqs(b,a,linspace(0,2*pi,100));

%% Exercise 2
Op = tan(0.4*pi/2);
Os = tan(0.5*pi/2);

T = 2;
Rs = 40;
Rp = 1;

[N, ws] = buttord(Op,Os,Rp,Rs,"s");
[bs,as] = butter(N,ws,"low","s");

[b,a] = bilinear(bs,as,1/T);
[B,W] = freqz(b,a,linspace(0,pi,512));

%% Exercise 3

wp = 0.4*pi;
ws = 0.5*pi;

Op = tan(wp/2);
Os = tan(ws/2);

Rp = 1;
Rs = 40;

[N, wc] = buttord(Op,Os,Rp,Rs,'s');
[bs,as] = butter(N,wc,'s');

T = 2;

[b,a] = bilinear(bs,as,1/T);

[B,W] = freqz(b,a,linspace(0,pi,512));

%% Exercise 4
wc = (pi/2)/pi;
Rp = 1;
Rs = 40;
N = 6;

[b1,a1] = butter(N,wc);
[b2,a2] = cheby1(N,Rp,wc);
[b3,a3] = cheby2(N,Rs,wc);
[b4,a4] = ellip(N,Rp,Rs,wc);

[B1,w1] = freqz(b1,a1,linspace(0,pi,512));
[B2,w2] = freqz(b2,a2,linspace(0,pi,512));
[B3,w3] = freqz(b3,a3,linspace(0,pi,512));
[B4,w4] = freqz(b4,a4,linspace(0,pi,512));

hold on
plot(w1/pi,20*log10(abs(B1)));
plot(w2/pi,20*log10(abs(B2)));
plot(w3/pi,20*log10(abs(B3)));
plot(w4/pi,20*log10(abs(B4)));
hold off
