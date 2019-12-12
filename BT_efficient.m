clear
clc
Nfft = 2048;
N = 1000;
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
noise_variance = 0.5:0.2:1.5;
no_trial = 50;
blackman = zeros(length(noise_variance),no_trial,fs/2);
blackman_mean = zeros(length(noise_variance),fs/2);
blackman_variance = zeros(length(noise_variance),fs/2);
blackman_bias = zeros(length(noise_variance),fs/2);
for step = 1:1:length(noise_variance)
    % Generate sample with white gaussian noise N(0,1)
    for trial = 1:1:no_trial
        w = sqrt(noise_variance(step))*randn(1,N);
        for n = 1:1:N
            x(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs)+w(n);
        end
        % Calculate correlation of signal
        Rx_bt = zeros(1,N);
        for m = (N/2):(N-1)
            for n = (N/2):(N-m-1+N/2)
                Rx_bt(m) = Rx_bt(m) + x(n+m-N/2)*x(n);
            end
            Rx_bt(m) = Rx_bt(m)/(N-m);
        end
        for m = 1:N/2
            Rx_bt(m) = Rx_bt(N-m);
        end
        % Generate Barlett window values
        M_bt = round(N/10);
        L = 2*M_bt+1;
        wbarlett0 = zeros(1,L);
        for m = 1:L
            if(m<=((L+1)/2))
                wbarlett0(m) = 2*m/(L-1);
            else
                wbarlett0(m) = 2-2*m/(L-1);
            end
        end
        % Calculating periodogram using Blackman-Turkey method with Barlett window
        P_bt = zeros(1,fs/2);
        for f = 1:fs/2
            for i = -M_bt:M_bt
                P_bt(f) = P_bt(f) + wbarlett0(i+(L+1)/2)*Rx_bt(N/2+i)*exp(-1i*2*pi*(f-1)*(i+N/2)*T);
            end
            P_bt(f) = abs(P_bt(f))*T;
        end
        P_bt = 10*log10(P_bt);
        for freq = 1:1:fs/2
            blackman(step,trial,freq) = P_bt(freq);
        end
    end
end

for step = 1:1:length(noise_variance)
    for trial = 1:1:no_trial
        for freq = 1:1:fs/2
            blackman_mean(step,freq) = blackman_mean(step,freq) + blackman(step,trial,freq);
        end
    end
    for freq = 1:1:fs/2
        blackman_mean(step,freq) = blackman_mean(step,freq)/no_trial;
    end
end

for step = 1:1:length(noise_variance)
    for trial = 1:1:no_trial
        for freq = 1:1:fs/2
            blackman_variance(step,freq) = (blackman_mean(step,freq) - blackman(step,trial,freq))^2 + blackman_variance(step,freq);
        end
    end
    for freq = 1:1:fs/2
        blackman_variance(step,freq) = blackman_variance(step,freq)/no_trial;
    end
end


% No noise
for n = 1:1:N
    x0(n) = a*cos(alpha*fs*2*pi*(n-1)/fs)+b*cos(beta*fs*2*pi*(n-1)/fs)+c*cos(muy*fs*2*pi*(n-1)/fs);
end

% Calculate correlation of signal
Rx0 = zeros(1,N);
for m = (N/2):(N-1)
    for n = (N/2):(N-m-1+N/2)
        Rx0(m) = Rx0(m) + x0(n+m-N/2)*x0(n);
    end
    Rx0(m) = Rx0(m)/(N-m);
end
for m = 1:N/2
    Rx0(m) = Rx0(N-m);
end

% Generate Barlett window values
M = round(N/10);
L = 2*M+1;
wbarlett0 = zeros(1,L);
for m = 1:L
    if(m<=((L+1)/2))
        wbarlett0(m) = 2*m/(L-1);
    else
        wbarlett0(m) = 2-2*m/(L-1);
    end
end

% Calculating periodogram using Blackman-Turkey method with Barlett window
Pbt0 = zeros(1,fs/2);
for f = 1:1:fs/2
    for i = -M:M
        Pbt0(f) = Pbt0(f) + wbarlett0(i+(L+1)/2)*Rx0(N/2+i)*exp(-1i*2*pi*(f-1)*(i+N/2)*T);
    end
    Pbt0(f) = abs(Pbt0(f))*T;
end
Pbt0 = 10*log10(Pbt0);

for step = 1:1:length(noise_variance)
    for freq = 1:1:fs/2
        blackman_bias(step,freq) = blackman_mean(step,freq)-Pbt0(freq);
    end
end

surf(blackman_mean)
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
surf(blackman_variance)
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
surf(blackman_bias)
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