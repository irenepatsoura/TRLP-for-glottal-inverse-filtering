[speech, fs] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav'); 
if size(speech, 2) > 1
    speech = mean(speech, 2);
end
speech = speech(:);

% Optional: same high-pass as in other experiments
fc_hp = 50;                         % Hz
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
speech = filtfilt(b_hp, a_hp, speech);

% -------------------------------------------------------------------------
% 2. Frame / LPC parameters
% -------------------------------------------------------------------------
frame_length = 0.3;    % seconds
frame_shift  = 0.02;   % seconds
p       = 20;          % LPC order
lambda1 = 1.0;         % diag regularization
lambda2 = 0.9;         % previous-coeff regularization
rho     = 0.99;        % leak factor for residual integration

frame_len_samples   = round(frame_length * fs);
frame_shift_samples = round(frame_shift * fs);
num_frames = floor((length(speech) - frame_len_samples) / frame_shift_samples) + 1;

% -------------------------------------------------------------------------
% 3. Buffers
% -------------------------------------------------------------------------
reconstructed_speech = zeros(length(speech), 1);   % synthesized speech
g_ola        = zeros(length(speech), 1);          % OLA glottal-like signal
window_sum   = zeros(length(speech), 1);          % OLA normalization

a_prev = zeros(p, 1);                             % previous LPC coeffs
win = hann(frame_len_samples, 'periodic');        % synthesis window

% -------------------------------------------------------------------------
% 4. Frame-by-frame TRLP + glottal estimation
% -------------------------------------------------------------------------
for k = 1:num_frames
    start_idx = (k - 1) * frame_shift_samples + 1;
    end_idx   = start_idx + frame_len_samples - 1;
    if end_idx > length(speech)
        break;
    end

    frame = speech(start_idx:end_idx);

    % --- Regularized autocorrelation LPC ---
    R = xcorr(frame, p, 'biased');
    R = R(p + 1:end);
    R_matrix = toeplitz(R(1:p));
    r_vector = R(2:p + 1);

    a_reg = lambda2 * a_prev;
    R_reg = R_matrix + lambda1 * eye(p);
    r_reg = r_vector + lambda1 * a_reg;

    a = R_reg \ r_reg;
    a_prev = a;

    % --- Inverse filtering: residual ---
    residual = filter([0; -a], 1, frame);

    % --- Simple leaky integration of residual -> glottal-like g_frame ---
    g_frame = filter(1, [1 -rho], residual);

    % --- Speech synthesis (reference only) ---
    synthesized_frame = filter(1, [1; -a], residual);

    % --- Overlap-add for speech ---
    reconstructed_speech(start_idx:end_idx) = ...
        reconstructed_speech(start_idx:end_idx) + synthesized_frame;

    % --- Overlap-add for glottal-like signal ---
    g_win = g_frame .* win;
    g_ola(start_idx:end_idx) = g_ola(start_idx:end_idx) + g_win;
    window_sum(start_idx:end_idx) = window_sum(start_idx:end_idx) + win;
end

% Normalize OLA glottal signal by window sum
window_sum(window_sum < 1e-3) = 1;
g_ola = g_ola ./ window_sum;

% Normalize reconstructed speech (for listening, optional)
reconstructed_speech = reconstructed_speech / max(abs(reconstructed_speech) + eps);
audiowrite('reconstructed_speech_trlp.wav', reconstructed_speech, fs);

% -------------------------------------------------------------------------
% 5. Simple speech comparison plot (optional)
% -------------------------------------------------------------------------
figure;
t_speech = (0:length(speech)-1)/fs;
subplot(2,1,1);
plot(t_speech, speech);
title('Original Speech Signal');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2);
plot(t_speech, reconstructed_speech);
title('Reconstructed Speech Signal (TRLP)');
xlabel('Time (s)'); ylabel('Amplitude');
set(gcf, 'Color', 'w');

% -------------------------------------------------------------------------
% 6. Load ground-truth glottal flow (ug)
% -------------------------------------------------------------------------
gt_filepath = 'data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav';
fprintf('Loading ground truth from: %s\n', gt_filepath);
[gt, fs_gt] = audioread(gt_filepath);
if size(gt, 2) > 1
    gt = mean(gt, 2);
end
gt = gt(:);

if fs_gt ~= fs
    fprintf('Resampling ground truth from %d Hz to %d Hz\n', fs_gt, fs);
    gt = resample(gt, fs, fs_gt);
end

% -------------------------------------------------------------------------
% 7. Align lengths, detrend, normalize
% -------------------------------------------------------------------------
min_length = min(length(gt), length(g_ola));
gt = gt(1:min_length);
g_ola = g_ola(1:min_length);

gt = detrend(gt, 'linear');
g_ola = detrend(g_ola, 'linear');

gt_norm = gt / (max(abs(gt)) + eps);
g_norm  = g_ola / (max(abs(g_ola)) + eps);

% -------------------------------------------------------------------------
% 8. Metrics
% -------------------------------------------------------------------------
correlation = corr(gt_norm, g_norm);
mse  = mean((gt_norm - g_norm).^2);
rmse = sqrt(mse);

fprintf('\n=== TRLP vs Ground Truth Metrics ===\n');
fprintf('Correlation: %.4f\n', correlation);
fprintf('RMSE:       %.4f\n', rmse);
fprintf('====================================\n\n');

% -------------------------------------------------------------------------
% 9. Plots (same style as QCP comparison)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 1600, 900]);
set(gcf, 'Color', 'w');

t = (0:min_length-1) / fs;

% 1) Full signal comparison
subplot(4, 2, [1, 2]);
plot(t, gt_norm, 'b', 'LineWidth', 1, 'DisplayName', 'Ground Truth');
hold on;
plot(t, g_norm, 'r', 'LineWidth', 1, 'DisplayName', 'TRLP Estimated');
xlabel('Time (s)');
ylabel('Normalized Amplitude');
title('Glottal Flow Comparison - TRLP');
legend('Location', 'best');
grid on;
hold off;

% 2) Zoom: beginning (0–50 ms)
subplot(4, 2, 3);
zoom_samples = round(0.05 * fs);  % first 50 ms
zoom_idx = 1:min(zoom_samples, min_length);
plot(t(zoom_idx), gt_norm(zoom_idx), 'b', 'LineWidth', 1.5);
hold on;
plot(t(zoom_idx), g_norm(zoom_idx), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoom: Beginning (0–50 ms)');
legend('Ground Truth', 'TRLP');
grid on;

% 3) Zoom: middle
subplot(4, 2, 4);
mid_start = round(min_length/2);
zoom_idx = mid_start:min(mid_start + zoom_samples, min_length);
plot(t(zoom_idx), gt_norm(zoom_idx), 'b', 'LineWidth', 1.5);
hold on;
plot(t(zoom_idx), g_norm(zoom_idx), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoom: Middle');
legend('Ground Truth', 'TRLP');
grid on;

% 4) Error signal
subplot(4, 2, [5, 6]);
error_signal = gt_norm - g_norm;
plot(t, error_signal, 'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Error');
title(sprintf('Estimation Error (RMSE: %.4f)', rmse));
grid on;

% 5) Scatter plot
subplot(4, 2, 7);
scatter(gt_norm, g_norm, 1, 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
plot([-1, 1], [-1, 1], 'r--', 'LineWidth', 2);  % ideal line
xlabel('Ground Truth');
ylabel('TRLP Estimated');
title(sprintf('Scatter Plot (Corr: %.4f)', correlation));
axis equal;
grid on;
xlim([-1, 1]);
ylim([-1, 1]);

% 6) Frequency-domain comparison (PSD)
subplot(4, 2, 8);
nfft = 2048;
[Pgt, f]  = pwelch(gt_norm, hann(512), 256, nfft, fs);
[Pest, ~] = pwelch(g_norm,  hann(512), 256, nfft, fs);

plot(f, 10*log10(Pgt), 'b', 'LineWidth', 1.5, 'DisplayName', 'Ground Truth');
hold on;
plot(f, 10*log10(Pest), 'r', 'LineWidth', 1.5, 'DisplayName', 'TRLP');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density');
legend('Location', 'best');
grid on;
xlim([0, 2000]);  % focus on lower frequencies