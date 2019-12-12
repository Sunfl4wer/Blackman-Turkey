clear;
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
variance = 0;
w = sqrt(variance)*randn(1,N);
% Generate sample with white gaussian noise N(0,1)
for n = 1:1:N
    x(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs)+w(n);
end

% Calculate correlation of signal
rx = zeros(1,N);
for m = (N/2):(N-1)
    for n = (N/2):(N-m-1+N/2)
        rx(m) = rx(m) + x(n+m-N/2)*x(n);
    end
    rx(m) = rx(m)/(N-m);
end
for m = 1:N/2
    rx(m) = rx(N-m);
end

% Generate Barlett window values
M = round(N/10);
L = 2*M+1;
wbarlett = zeros(1,L);
for m = 1:L
    if(m<=((L+1)/2))
        wbarlett(m) = 2*m/(L-1);
    else
        wbarlett(m) = 2-2*m/(L-1);
    end
end

% Calculating periodogram using Blackman-Turkey method with Barlett window
Pbt = zeros(1,fs/2);
for f = 1:fs/2
    for i = -M:M
        Pbt(f) = Pbt(f) + wbarlett(i+(L+1)/2)*rx(N/2+i)*exp(-1i*2*pi*(f-1)*(i+N/2)*T);
    end
    Pbt(f) = abs(Pbt(f))*T;
end
Pbt = 10*log10(Pbt);

% Plotting the signal and periodogram
subplot(1,2,1);
plot(t,x);
title_1 = strcat("Sample signal with f1 = ",num2str(f1),"Hz, f2 = ",num2str(f2),"Hz, f3 = ",num2str(f3),"Hz, with fs = ",num2str(fs),"Hz");
title_11 = strcat("number of samples = ",num2str(N),", noise = N(0,1)");
title({title_1,title_11});
subplot(1,2,2);
plot(0:1:fs/2-1,Pbt);
hold on
[pks,locs] = findpeaks(Pbt,0:1:fs/2-1);
findpeaks(Pbt,0:1:fs/2-1)
for i = 1:1:length(pks)
    text(locs(i)+.02,pks(i),num2str(locs(i)));
end
title("Mordified Periodogram - Blackman-Turkey");