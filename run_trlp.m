% Frame-by-frame TRLP Analysis
addpath("core","eval","framework","pipeline")

% 1. Load speech signal
[x, fs] = audioread('data/Glottal_signals_db/aa-pout-105Hz-8kHz.wav');

if size(x, 2) > 1
    x = mean(x, 2);
end
x = x(:);  % Ensure column vector

% --- Fix Polarity ---
% Glottal pulses have positive skewness. If negative, flip the signal.
if skewness(x) < 0
    x = -x;
    fprintf('Signal polarity inverted based on skewness check.\n');
end

% High-pass filtering
fc_hp = 50;
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
x = filtfilt(b_hp, a_hp, x);

% 2. Frame-by-frame analysis parameters
analysis_mode = 'fixed';
frame_length_ms = 50;
frame_shift_ms = 10;
frame_shift = round(frame_shift_ms * fs / 1000);
frame_length = round(frame_length_ms * fs / 1000);

% 3. Create frames
fprintf('Creating frames...\n');
[frames, frame_indices] = create_fixed_frames(x, frame_length, frame_shift);
num_frames = length(frames);
fprintf('Total frames: %d\n', num_frames);

% 4. TRLP options
p = 20;          % LPC order
rho = 0.99;      % Integrator leak factor
a_prev = zeros(p, 1); % Initialize previous coefficients

% 5. Process each frame
fprintf('Processing frames with TRLP...\n');
results = cell(num_frames, 1);

for i = 1:num_frames
    if mod(i, 10) == 0
        fprintf('  Frame %d/%d\n', i, num_frames);
    end
    
    frame_data = frames{i};
    
    % Apply TRLP to frame
    try
        % Calculate LPC coefficients
        a = mytrlp(frame_data, p, a_prev);
        
        % Update previous coefficients for next frame
        a_prev = a;
        
        % Inverse filtering to get residual
        % Note: filter([1; -a], 1, frame_data) corresponds to A(z) = 1 - sum(a_k z^-k)
        % But usually LPC returns a such that A(z) = 1 + sum(a_k z^-k) if using poly(a)
        % mytrlp returns 'a' from R \ r.
        % If x(n) = sum(a_k x(n-k)) + e(n), then e(n) = x(n) - sum(a_k x(n-k))
        % In filter terms: y = filter([1; -a], 1, x)
        residual = filter([1; -a], 1, frame_data);
        
        % Integration to get glottal flow
        g = filter(1, [1 -rho], residual);
        
        % Store results
        results{i}.g = g(:);
        results{i}.a = a;
        results{i}.frame_idx = frame_indices(i);
        results{i}.frame_length = length(frame_data);
    catch ME
        warning('Frame %d failed: %s', i, ME.message);
        results{i} = [];
    end
end

fprintf('Done processing!\n');

% 6. Compare with ground truth
gt_file = 'data/Glottal_signals_db/aa-ug-105Hz-8kHz.wav';
compare_with_ground_truth(x, fs, results, frame_indices, 'TRLP', gt_file);
