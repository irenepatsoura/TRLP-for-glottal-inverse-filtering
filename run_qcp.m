% 1. Load or create a speech signal
addpath("core","eval","framework","pipeline")
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

fc_hp = 50;  % Hz
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
x = filtfilt(b_hp, a_hp, x);

% 2. Frame-by-frame analysis parameters
analysis_mode = 'fixed';  % 'fixed' or 'adaptive'
frame_length_ms = 50;     % Frame length in ms (for fixed mode)
frame_shift_ms = 10;      % Frame shift in ms (50% overlap is common)
frame_shift = round(frame_shift_ms * fs / 1000);

% 3. Estimate pitch for adaptive mode
if strcmp(analysis_mode, 'adaptive')
    fprintf('Estimating pitch...\n');
    pitch_info = estimate_pitch(x, fs);
else
    frame_length = round(frame_length_ms * fs / 1000);
    pitch_info = struct('f0', [], 'periods', []);
end

% 3b. Choose a global f0 for QCP
global_f0 = 105;

% 4. Create frames
fprintf('Creating frames...\n');
if strcmp(analysis_mode, 'fixed')
    [frames, frame_indices] = create_fixed_frames(x, frame_length, frame_shift);
else
    [frames, frame_indices] = create_adaptive_frames(x, fs, pitch_info, frame_shift);
end

num_frames = length(frames);
fprintf('Total frames: %d\n', num_frames);



% 5. QCP options
options = struct();
options.dq = 0.4; % Reduced to ensure closed phase analysis
options.pq = 0.05; % Start window shortly after GCI
options.nramp = round(fs/8000*7);
options.causality = 0; % 0=noncausal integration (backward), 1=causal
options.f0 = global_f0;
options.remove_real_poles = 1; % Force removal of spectral tilt (source) from VT model

% 6. Process each frame
fprintf('Processing frames...\n');
results = cell(num_frames, 1);

for i = 1:num_frames
    if mod(i, 10) == 0
        fprintf('  Frame %d/%d\n', i, num_frames);
    end
    
    % Create signal object for current frame
    frame_signal = signal(frames{i}, fs);
    
    % Apply QCP to frame
    try
        [sg, Hvt, e_ar, Hg] = qcp(frame_signal, options);
        
        % Store results
        results{i}.sg = sg.s;
        results{i}.Hvt = Hvt;
        results{i}.e_ar = e_ar;
        results{i}.Hg = Hg;
        results{i}.frame_idx = frame_indices(i);
        results{i}.frame_length = length(frames{i});
    catch ME
        warning('Frame %d failed: %s', i, ME.message);
        results{i} = [];
    end
end

fprintf('Done processing!\n');

% 8. Compare with ground truth
gt_file = 'data/Glottal_signals_db/aa-ug-105Hz-8kHz.wav';  % Adjust path
compare_with_ground_truth(x, fs, results, frame_indices, 'QCP', gt_file);