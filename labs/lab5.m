%% Exercise 1
h1 = [1 2 3 -2 5 -2 3 2 1];
h2 = [1 2 3 -2 -2 3 2 1];
h3 = [1 2 3 -2 0 2 -3 -2 -1];
h4 = [1 2 3 -2 2 -3 -2 -1];

[H1,w1] = freqz(h1,1,'whole');
[H2,w2] = freqz(h2,1,'whole');
[H3,w3] = freqz(h3,1,'whole');
[H4,w4] = freqz(h4,1,'whole');

% hold on
% subplot(3,1,1)
% sgtitle('Magnitude and phase responses')
% plot(w1/pi, [abs(H1),abs(H2),abs(H3),abs(H4)]);
% xlabel('w/pi');
% ylabel('|H_i|');
% legend('H1','H2','H3','H4')
%
% subplot(3,1,2)
% plot(w1/pi, [angle(H1),angle(H2),angle(H3),angle(H4),]);
% xlabel('w/pi');
% ylabel('angle(H_i)');
% legend('H1','H2','H3','H4')
% hold off
%

%% Exercise 2
h = [-4 1 -1 -2 5 6 6 5 -2 -1 1 -4];

[A, w, phi] = zerophase(h,1,"whole");

% subplot(311)
% sgtitle('Magnitude and phaseresponse of H')
% xlabel('w/pi')
% ylabel('|H(jw)|')
% plot(w/pi,A)
% yline(0)
% subplot(312)
% xlabel('w/pi')
% ylabel('Angle(H)')
% plot(w/pi,phi)
% subplot(313)
% zplane(h)

%% Exercise 3

h0 = [2];
h1 = [1.5];
h2 = [-3.2];
h3 = [-5.2];
h4 = [6.4];

h = [h0 h1 h2 h3 h4];

% Type I
h_I = [h fliplr(h(1:end-1))];

% Type II
h_II = [h fliplr(h)];

% Type III
h_III = [h 0 -fliplr(h)];

% Type IV
h_IV = [h -fliplr(h)];

% sgtitle('Zeros of filters')
% subplot(221)
% zplane(h_I)
% title('h_I')

% subplot(222)
% zplane(h_II)
% title('h_{II}')

% subplot(223)
% zplane(h_III)
% title('h_{III}')

% subplot(224)
% zplane(h_IV)
% title('h_{IV}')

%% Exercise 4

N = 50;
n = -N/2:N/2;

wc = (pi/2)/pi;
h_lp = wc*sinc(wc*n);
[H_lp,w] = freqz(h_lp,1,linspace(-pi,pi,512));

% sgtitle('Truncated lowpassfilter')
% subplot(311)
% plot(w/pi,abs(H_lp))
% title('Magnitude |H_{lp}|')
% xlabel('w/pi')
% ylabel('Magnitude of h_{lp}')
% subplot(312)
% plot(w/pi,angle(H_lp))
% title('Phase of H_{lp}')
% xlabel('w/pi')
% ylabel('Angle(H_{lp})')
% subplot(313)
% zplane(h_lp)

%% Exercise 5
N = 51;


h_rect = ones(1,N);
h_hamming = hamming(N);
h_blackman = blackman(N);
wtest = linspace(-pi,pi,1001);

[H_rect,wrect] = zerophase(h_rect,1,w);
[H_hamming, whamming] = zerophase(h_hamming,1,w);
[H_blackman, wblackman] = zerophase(h_blackman,1,w);


% sgtitle('Amplitudes of windows')
% subplot(211)
% title('Amplitudes')
% hold on
% xlabel('w/pi')
% ylabel('|H_i(jw)|')
% plot(w/pi,(H_rect))
% plot(w/pi,(H_hamming))
% plot(w/pi,(H_blackman))
% legend('Rectangular', 'Hamming', 'Blackman')
% hold off
%
% subplot(212)
% title('Amplitudes in dB')
% hold on
% xlabel('w/pi')
% ylabel('|H_i(jw)| in dB')
% plot(w/pi,20*log10(abs(H_rect)/max(abs(H_rect))))
% plot(w/pi,20*log10(abs(H_hamming)/max(abs(H_hamming))))
% plot(w/pi,20*log10(abs(H_blackman)/max(abs(H_blackman))))
% legend('Rectangular', 'Hamming', 'Blackman')
% hold off

%% Exercise 6

wc = (pi/2)/pi;
w = linspace(-pi,pi,1024);
N = 50;
n = -N/2:N/2;

h_rect = ones(1,N+1);
h_hamming = hamming(N+1);
h_blackman = blackman(N+1);

h_lp = wc*sinc(wc*n);

% Apply window
x1 = h_rect.*h_lp;
x2 = h_hamming.*h_lp;
x3 = h_blackman.*h_lp;

[H_lp,wh] = zerophase(h_lp,1,w);
[X1,w1] = zerophase(x1,1,w);
[X2,w2] = zerophase(x2,1,w);
[X3,w3] = zerophase(x3,1,w);

sgtitle('Amplitude responses')
subplot(211)
title('Ideal lowpass without specific windows')
plot(wh/pi,H_lp)
xlabel('wh/pi')
ylabel('|H_{lp}(jw)|')
subplot(212)
title('Lowpass windowed with specific windows')
hold on
plot(w1/pi,(X1))
plot(w2/pi,(X2))
plot(w3/pi,(X3))
hold off
xlabel('w/pi')
ylabel('|X_i(jw)|')