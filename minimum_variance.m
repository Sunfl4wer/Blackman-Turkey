clear,clc;
Nfft = 2048;
N = 100;
a = 5;
b = 2;
c = 5;
alpha = 0.37;
beta = 0.25;
muy = 0.4;
fs = 2000;
T = 1/fs;
t = 0:1:N-1;
t = t/fs;
f1 = alpha*fs;
f2 = beta*fs;
f3 = muy*fs;
variance = 1;
stddv= sqrt(variance);
w = stddv*randn(1,N);
% Generate sample with white gaussian noise N(0,1)
for n = 1:1:N
    sig(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs)+w(n);
end
for n = 1:1:N
    xm(n) = sig(n);
end
for n = N+1:1:Nfft
    xm(n) = 0;
end
xm=xm-mean(xm);    % remove mean to prepare for covariance estimation
xc=xcorr(xm,Nfft-1,'biased'); % compute covariance sequence
Rx=toeplitz(xc(Nfft:end));
v = zeros(Nfft,fs);
for j = 1:1:fs
    for q = 1:1:Nfft
        v(q,j) = exp(-1i*2*pi*(j-1)*(q-1)*T);
    end
end
Px = (v')*(inv(Rx))*(v);
P = zeros(1,fs/2);
for q = 1:1:fs/2
        P(q) = Px(q,q);
end
for q = 1:1:fs/2
        P(q) = abs(P(q));
        P(q) = Nfft/P(q);
        P(q) = 10*log10(P(q));
end

% Plotting the signal and periodogram
subplot(1,2,1);
plot(t,sig),grid on;
title_1 = strcat("Sample signal with f1 = ",num2str(f1),"Hz, f2 = ",num2str(f2),"Hz, f3 = ",num2str(f3),"Hz, with fs = ",num2str(fs),"Hz");
title_11 = strcat("number of samples = ",num2str(N),", noise = N(0,1)");
title({title_1,title_11});
subplot(1,2,2);
plot(0:1:fs/2-1,P),grid on;
hold on
[pks,locs] = findpeaks(P,0:1:fs/2-1);
findpeaks(P,0:1:fs/2-1)
for i = 1:1:length(pks)
    text(locs(i)+.02,pks(i),num2str(locs(i)));
end
title("Minimum variance");
hold off