% ==== FULL-SIGNAL IAIF TEST ====

% 1. Load and preprocess speech (same as for frames)
[x, fs] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');
x = mean(x, 2);
x = x(:);

% High-pass pre-emphasis / DC removal (IAIF block 1)
fc_hp = 50;
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
x = filtfilt(b_hp, a_hp, x);

s_full = signal(x, fs);

% 2. IAIF options (same as frame script)
options = struct();
options.p = 20;
options.g = 4;
options.f0 = 105;      % known F0
options.causality = 0; % noncausal
% options.rho = 0.99;  % optional

% 3. Run IAIF on the FULL signal
[g_full, Hvt_full, e_ar_full, Hg_full] = iaif(s_full, options);

if isa(g_full, 'signal')
    g_est = g_full.s(:);
else
    g_est = g_full(:);
end

% 4. Load ground-truth glottal flow (ug)
[ug, fs_ug] = audioread( ...
    'data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');
if fs_ug ~= fs
    error('Sampling rate mismatch between speech and ug.');
end
ug = ug(:);

% 5. Make same length
N = min(length(g_est), length(ug));
g_est = g_est(1:N);
ug = ug(1:N);

% 6. Detrend both
g_est = detrend(g_est, 'linear');
ug = detrend(ug, 'linear');

% 7. Corr before any lag correction
c0 = corr(g_est, ug);
fprintf('Corr(IAIF g, ug) before lag correction: %f\n', c0);

% 8. Optional: estimate constant lag (like we did for QCP)
max_lag_s = 0.010;         % +/- 10 ms
max_lag = round(max_lag_s * fs);

[cc, lags] = xcorr(g_est, ug, max_lag, 'coeff');
[~, idx_max] = max(cc);
best_lag = lags(idx_max);
fprintf('Best lag (g vs ug): %d samples (%.4f s)\n', ...
        best_lag, best_lag/fs);

if best_lag > 0
    % g_est is ahead
    g_sync = g_est(1+best_lag:end);
    ug_sync = ug(1:end-best_lag);
elseif best_lag < 0
    shift = -best_lag;
    g_sync = g_est(1:end-shift);
    ug_sync = ug(1+shift:end);
else
    g_sync = g_est;
    ug_sync = ug;
end

Ns = min(length(g_sync), length(ug_sync));
g_sync = g_sync(1:Ns);
ug_sync = ug_sync(1:Ns);

c1 = corr(g_sync, ug_sync);
fprintf('Corr(IAIF g, ug) after lag correction:  %f\n', c1);

% 9. Plots
t = (0:Ns-1)/fs;

figure;
plot(t, ug_sync / max(abs(ug_sync)), 'b'); hold on;
plot(t, g_sync / max(abs(g_sync)), 'r--');
legend('Ground truth ug', 'IAIF g');
xlabel('Time (s)');
ylabel('Normalized amplitude');
title('Full-signal IAIF vs Ground truth ug (aligned)');
grid on;

% Zoom 0.20–0.23 s
t_start = 0.20; t_end = 0.23;
idx = round(t_start*fs):round(t_end*fs);
idx = idx(idx >= 1 & idx <= Ns);

figure;
plot(t(idx), ug_sync(idx)/max(abs(ug_sync(idx))), 'b'); hold on;
plot(t(idx), g_sync(idx)/max(abs(g_sync(idx))), 'r--');
legend('Ground truth ug', 'IAIF g');
xlabel('Time (s)');
ylabel('Normalized amplitude');
title(sprintf('Zoomed IAIF vs ug (%.3f–%.3f s)', t_start, t_end));
grid on;