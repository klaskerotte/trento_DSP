%% Task 1
% Plot x[n] and the magnitude of its DFT |X(jw)|
function [delta, n] = impseq(n0,nstart,nend)
n = nstart:nend;
delta = [(n==n0)];
end

N = 200;
n = 0:N-1;

x = zeros(1,N);
for r = 1:25
    x = x + 1/25.*impseq(8*r-1,0,N-1);
end
% plot(n,x)

%%  Task 2
% Plot the magnitude of the DFT
[X,w] = freqz(x,1,1024);
% plot(w/pi,abs(X));

%% Task 3
% Create a filter that attenuates the three first peaks by at least 65dB while not attenuating the other peaks more than 3dB

M = 500;
w_cutoff = 0.65;
w_lastpeak = 0.5;
w_no_touch = 0.75;
[b,a] = fir1(M,w_cutoff,"high");
[B,w] = freqz(b,a,1024);


% Using Parks-McLellan:
% MP = 200;
% transitionbands = [0.6 0.68];
% fp = [0 transitionbands 1];
% amps = [0 10^(-65/20) 1 1];
% bp = firpm(MP,fp,amps);
% [BP,w] = freqz(bp,1,1024);
%
% subplot(2,1,1)
% plot(w/pi,20*log10(abs(BP)))
% yline([-65,0])
% subplot(2,1,2)
% plot(w/pi,angle(BP))





% subplot(2,1,1)
% plot(w/pi,abs(X));
% subplot(2,1,2)
% plot(w/pi,20*log10(abs(B)))
% xline(w_cutoff)
% yline(-65)

%% Task 4
% Plot the outputvalues y[n] when the input is x[n]

yb = filter(b,a,x);


%% Task 5
d = 0.038;
[order, wn] = buttord(w_cutoff-d,w_cutoff+d,3,65)
[bfb,bfa] = butter(order,wn,"high");

[BF,w] = freqz(bfb,bfa,1024);
% plot(w/pi,20*log10(abs(BF)))
% xline(w_cutoff)
% xline(w_lastpeak)
% xline(w_no_touch)
% yline(-65)

%% Task 6

Y = BF.*X;
ybf = filter(b,a,x);

% subplot(2,1,1)
% plot(w/pi,abs(Y));
% xline(w_no_touch)
% xline(w_lastpeak)
% xline(w_cutoff)
% yline(-65)
% yline(-3)



