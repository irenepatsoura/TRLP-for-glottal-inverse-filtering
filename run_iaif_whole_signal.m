% Whole Signal IAIF Analysis (No Frames)
% For comparison with frame-by-frame approach

% 1. Load speech signal
[x, fs] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');

if size(x, 2) > 1
    x = mean(x, 2);
end
x = x(:);  % Ensure column vector

fprintf('Processing whole signal with IAIF...\n');
fprintf('Signal length: %.3f seconds (%d samples)\n', length(x)/fs, length(x));

% 2. IAIF options
options = struct();
options.p = 20;  % Vocal tract order
options.g = 4;   % Glottal source order

% 3. Apply IAIF to entire signal
[g, Hvt, e_ar, Hg] = iaif(x, fs, options);

% Handle signal object if needed
if isa(g, 'signal')
    g = g.s;
end

fprintf('Done processing!\n');

% 4. Visualize whole-signal results
visualize_whole_signal_iaif(x, fs, g, Hvt, e_ar, Hg);

% 5. Compare with ground truth
gt_file = 'data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav';
compare_whole_signal_with_gt(x, fs, g, gt_file);

%% Visualization function for whole signal
function visualize_whole_signal_iaif(x, fs, g, Hvt, e_ar, Hg)
    figure('Position', [50, 50, 1600, 900]);
    
    % Time vectors
    t_x = (0:length(x)-1) / fs;
    t_g = (0:length(g)-1) / fs;
    t_ar = (0:length(e_ar)-1) / fs;
    
    % 1. Original signal
    subplot(4, 2, [1, 2]);
    plot(t_x, x, 'Color', [0.5, 0.5, 0.5]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original Speech Signal');
    grid on;
    
    % 2. Estimated glottal flow - full
    subplot(4, 2, [3, 4]);
    plot(t_g, g, 'b', 'LineWidth', 0.8);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Estimated Glottal Flow (Whole Signal)');
    grid on;
    
    % 3. Glottal flow - zoomed beginning
    subplot(4, 2, 5);
    zoom_samples = round(0.05 * fs);  % 50ms
    zoom_idx = 1:min(zoom_samples, length(g));
    plot(t_g(zoom_idx), g(zoom_idx), 'b', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Glottal Flow - Zoom: Beginning (0-50ms)');
    grid on;
    
    % 4. Glottal flow - zoomed middle
    subplot(4, 2, 6);
    mid_start = round(length(g)/2);
    zoom_idx = mid_start:min(mid_start + zoom_samples, length(g));
    plot(t_g(zoom_idx), g(zoom_idx), 'b', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Glottal Flow - Zoom: Middle');
    grid on;
    
    % 5. AR residual
    subplot(4, 2, 7);
    plot(t_ar, e_ar, 'Color', [0.8, 0.4, 0]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('AR Residual (Glottal Excitation)');
    grid on;
    
    % 6. Frequency responses
    subplot(4, 2, 8);
    [H_vt, f] = freqz(1, Hvt, 2048, fs);
    [H_g, ~] = freqz(1, Hg, 2048, fs);
    
    plot(f, 20*log10(abs(H_vt)), 'b', 'LineWidth', 1.5, 'DisplayName', 'Vocal Tract');
    hold on;
    plot(f, 20*log10(abs(H_g)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Glottal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('Filter Frequency Responses');
    legend('Location', 'best');
    grid on;
    xlim([0, 4000]);
    hold off;
    
    set(gcf, 'Color', 'w');
end

%% Comparison with ground truth for whole signal
function compare_whole_signal_with_gt(x, fs, g_estimated, gt_filepath)
    % Load ground truth
    fprintf('Loading ground truth from: %s\n', gt_filepath);
    [gt, fs_gt] = audioread(gt_filepath);
    
    if size(gt, 2) > 1
        gt = mean(gt, 2);
    end
    gt = gt(:);
    
    % Ensure column vector
    g_estimated = g_estimated(:);
    
    % Resample if needed
    if fs_gt ~= fs
        fprintf('Resampling ground truth from %d Hz to %d Hz\n', fs_gt, fs);
        gt = resample(gt, fs, fs_gt);
    end
    
    fprintf('Ground truth length: %d samples\n', length(gt));
    fprintf('Estimated length: %d samples\n', length(g_estimated));
    
    % Align signals - use minimum length
    min_length = min(length(gt), length(g_estimated));
    gt = gt(1:min_length);
    g_estimated = g_estimated(1:min_length);
    
    fprintf('Using %d samples for comparison\n', min_length);
    
    % Normalize both signals
    gt_norm = gt / max(abs(gt) + eps);
    g_norm = g_estimated / max(abs(g_estimated) + eps);
    
    % Ensure same size (sometimes there are dimension mismatches)
    gt_norm = gt_norm(:);
    g_norm = g_norm(:);
    
    % Double check sizes
    if length(gt_norm) ~= length(g_norm)
        error('Mismatch after normalization: GT=%d, EST=%d', length(gt_norm), length(g_norm));
    end
    
    % Calculate metrics
    correlation = corr(gt_norm, g_norm);
    mse = mean((gt_norm - g_norm).^2);
    rmse = sqrt(mse);
    
    fprintf('\n=== WHOLE SIGNAL Comparison Metrics ===\n');
    fprintf('Correlation: %.4f\n', correlation);
    fprintf('RMSE: %.4f\n', rmse);
    fprintf('=======================================\n\n');
    
    % Create comparison figure
    figure('Position', [100, 100, 1600, 900]);
    
    t = (0:min_length-1) / fs;
    
    % 1. Full signal comparison
    subplot(4, 2, [1, 2]);
    plot(t, gt_norm, 'b', 'LineWidth', 1, 'DisplayName', 'Ground Truth');
    hold on;
    plot(t, g_norm, 'r', 'LineWidth', 1, 'DisplayName', 'IAIF Estimated (Whole)');
    xlabel('Time (s)');
    ylabel('Normalized Amplitude');
    title('Glottal Flow Comparison - IAIF (Whole Signal)');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % 2. Zoomed view - beginning
    subplot(4, 2, 3);
    zoom_samples = round(0.05 * fs);
    zoom_idx = 1:min(zoom_samples, min_length);
    plot(t(zoom_idx), gt_norm(zoom_idx), 'b', 'LineWidth', 1.5);
    hold on;
    plot(t(zoom_idx), g_norm(zoom_idx), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Zoom: Beginning (0-50ms)');
    legend('Ground Truth', 'IAIF (Whole)');
    grid on;
    
    % 3. Zoomed view - middle
    subplot(4, 2, 4);
    mid_start = round(min_length/2);
    zoom_idx = mid_start:min(mid_start + zoom_samples, min_length);
    plot(t(zoom_idx), gt_norm(zoom_idx), 'b', 'LineWidth', 1.5);
    hold on;
    plot(t(zoom_idx), g_norm(zoom_idx), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Zoom: Middle');
    legend('Ground Truth', 'IAIF (Whole)');
    grid on;
    
    % 4. Error signal
    subplot(4, 2, [5, 6]);
    error_signal = gt_norm - g_norm;
    plot(t, error_signal, 'k', 'LineWidth', 0.5);
    xlabel('Time (s)');
    ylabel('Error');
    title(sprintf('Estimation Error (RMSE: %.4f)', rmse));
    grid on;
    
    % 5. Scatter plot
    subplot(4, 2, 7);
    scatter(gt_norm, g_norm, 1, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    plot([-1, 1], [-1, 1], 'r--', 'LineWidth', 2);
    xlabel('Ground Truth');
    ylabel('IAIF Estimated (Whole)');
    title(sprintf('Scatter Plot (Corr: %.4f)', correlation));
    axis equal;
    grid on;
    xlim([-1, 1]);
    ylim([-1, 1]);
    
    % 6. Frequency domain comparison
    subplot(4, 2, 8);
    nfft = 2048;
    [Pgt, f] = pwelch(gt_norm, hann(512), 256, nfft, fs);
    [Pest, ~] = pwelch(g_norm, hann(512), 256, nfft, fs);
    
    plot(f, 10*log10(Pgt), 'b', 'LineWidth', 1.5, 'DisplayName', 'Ground Truth');
    hold on;
    plot(f, 10*log10(Pest), 'r', 'LineWidth', 1.5, 'DisplayName', 'IAIF (Whole)');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density');
    legend('Location', 'best');
    grid on;
    xlim([0, 2000]);
    
    set(gcf, 'Color', 'w');
end