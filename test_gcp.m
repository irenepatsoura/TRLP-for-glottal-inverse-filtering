% Preprocess full speech (include your high?pass here)
[x, fs] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');
x = mean(x, 2);

% High-pass
fc_hp = 50;
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
x = filtfilt(b_hp, a_hp, x);

s_full = signal(x, fs);

options = struct();
options.dq = 0.7;
options.pq = 0.1;
options.nramp = round(fs/8000*7);
options.causality = 0;
options.f0 = 105;     % or from global pitch

[sg_full, Hvt_full, e_ar_full, Hg_full] = qcp(s_full, options);
g_est_full = sg_full.s(:);   % column  % make sure it is a column vector
dg_full = dg_full(:);

%% Load ground-truth glottal flow
% IMPORTANT: use the *-ug-*.wav file here, not the speech signal
[ug, fs_ug] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');
if fs_ug ~= fs
    error('Sampling rate mismatch between speech and ground-truth ug.');
end
ug = ug(:);

%% Make same length
N = min([length(g_est_full), length(dg_full), length(ug)]);
g_est_full = g_est_full(1:N);
dg_full    = dg_full(1:N);
ug = ug(1:N);

%% Detrend both in the same way
g_est_full = detrend(g_est_full, 'linear');
ug = detrend(ug, 'linear');

%% Correlation with ug vs derivative of ug
ug_diff = [0; diff(ug)];   % same length N

c1 = corr(g_est_full, ug);       % with glottal flow
c2 = corr(g_est_full, ug_diff);  % with derivative

fprintf('Correlation with ug:     %f\n', c1);
fprintf('Correlation with d(ug):  %f\n', c2);

% Optional: fix sign so we can at least compare shapes
if c1 < 0
    g_est_full = -g_est_full;
    c1 = -c1;
    c2 = -c2;
    fprintf('Sign flipped. New corr with ug: %f\n', c1);
end

%% Full-length normalized plot
t = (0:N-1) / fs;

figure;
plot(t, g_est_full / max(abs(g_est_full)), 'r--'); hold on;
plot(t, ug / max(abs(ug)), 'b');
legend('QCP estimate', 'Ground truth ug');
xlabel('Time (s)');
ylabel('Normalized amplitude');
title('Full-length comparison: QCP vs Ground truth');
grid on;

%% Zoomed plot over a few periods
t_start = 0.200;   % seconds
t_end   = 0.230;   % seconds
idx = round(t_start*fs) : round(t_end*fs);
idx(idx < 1 | idx > N) = [];

figure;
plot(t(idx), ug(idx) / max(abs(ug(idx))), 'b'); hold on;
plot(t(idx), g_est_full(idx) / max(abs(g_est_full(idx))), 'r--');
legend('Ground truth ug', 'QCP estimate');
xlabel('Time (s)');
ylabel('Normalized amplitude');
title(sprintf('Zoomed comparison (%.3f–%.3f s)', t_start, t_end));
grid on;