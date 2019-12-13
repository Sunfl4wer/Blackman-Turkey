clear
clc
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
f = 0:1:fs/2-1;
noise_variance = 0.5:0.2:1.5;
no_trial = 50;
mv = zeros(length(noise_variance),no_trial,fs/2);
mv_mean = zeros(length(noise_variance),fs/2);
mv_variance = zeros(length(noise_variance),fs/2);
mv_bias = zeros(length(noise_variance),fs/2);
sig = zeros(1,N);
xm = zeros(1,Nfft);
v = zeros(Nfft,fs);
for j = 1:1:fs
    for q = 1:1:Nfft
        v(q,j) = exp(-1i*2*pi*(j-1)*(q-1)*T);
    end
end
for step = 1:1:length(noise_variance)
    % Generate sample with white gaussian noise N(0,1)
    for trial = 1:1:no_trial
        w = sqrt(noise_variance(step))*randn(1,N);
        for n = 1:1:N
            xm(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs)+w(n);;
        end
        xm=xm-mean(xm);    % remove mean to prepare for covariance estimation
        xc=xcorr(xm,Nfft-1,'biased'); % compute covariance sequence
        Rx=toeplitz(xc(Nfft:end));
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
        for freq = 1:1:fs/2
            mv(step,trial,freq) = P(freq);
        end
    end
end

for step = 1:1:length(noise_variance)
    for trial = 1:1:no_trial
        for freq = 1:1:fs/2
            mv_mean(step,freq) = mv_mean(step,freq) + mv(step,trial,freq);
        end
    end
    for freq = 1:1:fs/2
        mv_mean(step,freq) = mv_mean(step,freq)/no_trial;
    end
end

for step = 1:1:length(noise_variance)
    for trial = 1:1:no_trial
        for freq = 1:1:fs/2
            mv_variance(step,freq) = (mv_mean(step,freq) - mv(step,trial,freq))^2 + mv_variance(step,freq);
        end
    end
    for freq = 1:1:fs/2
        mv_variance(step,freq) = mv_variance(step,freq)/no_trial;
    end
end


% No noise
for n = 1:1:N
    sig0(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs)+w(n);
end
for n = 1:1:N
    xm0(n) = sig0(n);
end
for n = N+1:1:Nfft
    xm0(n) = 0;
end
xm0=xm0-mean(xm0);    % remove mean to prepare for covariance estimation
xc0=xcorr(xm0,Nfft-1,'biased'); % compute covariance sequence
Rx0=toeplitz(xc0(Nfft:end));
Px0 = (v')*(inv(Rx0))*(v);
P0 = zeros(1,fs/2);
for q = 1:1:fs/2
        P0(q) = Px0(q,q);
end
for q = 1:1:fs/2
        P0(q) = abs(P0(q));
        P0(q) = Nfft/P0(q);
        P0(q) = 10*log10(P0(q));
end

for step = 1:1:length(noise_variance)
    mv_bias(step,freq) = mv_mean(step,freq)-P0(freq);
end


surf(mv_mean)
colormap(jet)
xlabel("Frequency(Hz)");
xticks(0:50:fs/2-1);
ylabel("Noise variance");
yticks(1:1:6);
yticklabels(noise_variance);
zlabel("PSD(dB)");
hold on
title_1 = strcat("Mean value of PSD with fs = ",num2str(fs),", AWGN with variance varies from 0.5 to 1.5, number of trial = ",num2str(no_trial));
title_2 = strcat("Signal has 3 components, f1 = ",num2str(f1),", f2 = ",num2str(f2),", f3 = ",num2str(f3),", number of samples is ",num2str(N));
title({title_1,title_2});
hold off

figure
surf(mv_variance)
colormap(jet)
xlabel("Frequency(Hz)");
xticks(0:50:fs/2-1);
ylabel("Noise variance");
yticks(1:1:6);
yticklabels(noise_variance);
zlabel("Variance(dB)");
hold on
title_1 = strcat("Variance value of PSD with fs = ",num2str(fs),", AWGN with variance varies from 0.5 to 1.5, number of trial = ",num2str(no_trial));
title_2 = strcat("Signal has 3 components, f1 = ",num2str(f1),", f2 = ",num2str(f2),", f3 = ",num2str(f3),", number of samples is ",num2str(N));
title({title_1,title_2});
hold off

figure
surf(mv_bias)
colormap(jet)
xlabel("Frequency(Hz)");
xticks(0:50:fs/2-1);
ylabel("Noise variance");
yticks(1:1:6);
yticklabels(noise_variance);
zlabel("bias(dB)");
hold on
title_1 = strcat("Bias value of PSD with fs = ",num2str(fs),", AWGN with variance varies from 0.5 to 1.5, number of trial = ",num2str(no_trial));
title_2 = strcat("Signal has 3 components, f1 = ",num2str(f1),", f2 = ",num2str(f2),", f3 = ",num2str(f3),", number of samples is ",num2str(N));
title({title_1,title_2});
hold off
