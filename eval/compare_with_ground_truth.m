function compare_with_ground_truth(x, fs, results, frame_indices, method_name, gt_filepath)
    % Simplified comparison: ground truth vs estimated (one above the other + zoom + PSD)
    
    % Determine field name
    if strcmp(method_name, 'QCP')
        field_name = 'sg';
    elseif strcmp(method_name, 'IAIF')
        field_name = 'g';
    else
        error('Unknown method. Use QCP or IAIF');
    end
    
    % Load ground truth
    fprintf('Loading ground truth from: %s\n', gt_filepath);
    [gt, fs_gt] = audioread(gt_filepath);
    if size(gt, 2) > 1, gt = mean(gt, 2); end
    gt = gt(:);
    
    if fs_gt ~= fs
        fprintf('Resampling ground truth from %d Hz to %d Hz\n', fs_gt, fs);
        gt = resample(gt, fs, fs_gt);
    end
    
    % Reconstruct estimated
    fprintf('Reconstructing estimated glottal flow...\n');
    estimated = reconstruct_signal(results, frame_indices, length(x), field_name, fs);
    estimated = detrend(estimated, 'linear');
    
    % Align
    min_length = min(length(gt), length(estimated));
    gt = gt(1:min_length);
    estimated = estimated(1:min_length);
    
    % Normalize
    gt_norm = gt / max(abs(gt));
    estimated_norm = estimated / max(abs(estimated));
    
    % Metrics
    correlation = corr(gt_norm, estimated_norm);
    rmse = sqrt(mean((gt_norm - estimated_norm).^2));
    
    fprintf('\n=== Comparison Metrics (%s) ===\n', method_name);
    fprintf('Correlation: %.4f\n', correlation);
    fprintf('RMSE:       %.4f\n', rmse);
    fprintf('================================\n\n');
    
    % === SIMPLIFIED PLOTS ===
    figure('Position', [100, 100, 1400, 600]);
    set(gcf, 'Color', 'w');
    
    t = (0:min_length-1) / fs;
    
    % 1) Ground truth (top)
    subplot(3, 2, [1, 2]);
    plot(t, gt_norm, 'b', 'LineWidth', 0.8);
    ylabel('Amplitude');
    title('Ground Truth Glottal Flow');
    grid on; xlim([0 t(end)]);
    
    % 2) Estimated (below ground truth)
    subplot(3, 2, [3, 4]);
    plot(t, estimated_norm, 'r', 'LineWidth', 0.8);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('%s Estimated Glottal Flow (Corr: %.3f, RMSE: %.3f)', ...
                  method_name, correlation, rmse));
    grid on; xlim([0 t(end)]);
    
    % 3) Zoom: 3-4 periods around middle
    subplot(3, 2, 5);
    % Pick a window of ~40 ms in the middle
    mid_t = t(end)/2;
    zoom_start_t = mid_t - 0.02;  % 20 ms before middle
    zoom_end_t   = mid_t + 0.02;  % 20 ms after
    zoom_idx = find(t >= zoom_start_t & t <= zoom_end_t);
    
    plot(t(zoom_idx), gt_norm(zoom_idx), 'b', 'LineWidth', 1.5, ...
         'DisplayName', 'Ground Truth');
    hold on;
    plot(t(zoom_idx), estimated_norm(zoom_idx), 'r--', 'LineWidth', 1.5, ...
         'DisplayName', method_name);
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Zoomed Comparison (3-4 periods)');
    legend('Location', 'best'); grid on;
    
    % 4) Power Spectral Density
    subplot(3, 2, 6);
    nfft = 2048;
    [Pgt, f]  = pwelch(gt_norm, hann(512), 256, nfft, fs);
    [Pest, ~] = pwelch(estimated_norm, hann(512), 256, nfft, fs);
    
    plot(f, 10*log10(Pgt), 'b', 'LineWidth', 1.5, 'DisplayName', 'Ground Truth');
    hold on;
    plot(f, 10*log10(Pest), 'r', 'LineWidth', 1.5, 'DisplayName', method_name);
    xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density');
    legend('Location', 'best'); grid on;
    xlim([0, 2000]);
end