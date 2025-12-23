% Frame-by-frame IAIF Analysis
% 1. Load speech signal
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

%%% NEW: High-pass filtering (IAIF block 1 must be done in advance)
fc_hp = 50;                          % e.g. 40�70 Hz
[b_hp, a_hp] = butter(2, fc_hp/(fs/2), 'high');
x = filtfilt(b_hp, a_hp, x);
%%% END NEW

% 2. Frame-by-frame analysis parameters
analysis_mode = 'fixed';  % 'fixed' or 'adaptive'
frame_length_ms = 50;     % Frame length in ms (for fixed mode)
frame_shift_ms = 10;      % Frame shift in ms
frame_shift = round(frame_shift_ms * fs / 1000);

% 3. Estimate pitch (only needed if you want adaptive f0 per frame)
if strcmp(analysis_mode, 'adaptive')
    fprintf('Estimating pitch (for adaptive f0)...\n');
    pitch_info = estimate_pitch(x, fs);  % must return .f0 per frame or sample
else
    frame_length = round(frame_length_ms * fs / 1000);
    pitch_info = struct('f0', [], 'periods', []);
end

% 3b. Choose a global f0 for IAIF (to avoid find_f0 on each 30-ms frame)
%%% NEW: choose global f0 once (here from file name, 105 Hz)
global_f0 = 105;  % or compute once from estimate_pitch(x, fs)
%%% END NEW

% 4. Create frames (reuses QCP functions)
fprintf('Creating frames...\n');
if strcmp(analysis_mode, 'fixed')
    [frames, frame_indices] = create_fixed_frames(x, frame_length, frame_shift);
else
    [frames, frame_indices] = create_adaptive_frames(x, fs, pitch_info, frame_shift);
end

num_frames = length(frames);
fprintf('Total frames: %d\n', num_frames);

% 5. IAIF options
options = struct();
options.p =  fs/1000 + 2;     % Vocal tract order
options.g = 4;      % Glottal source order

%%% NEW: pass f0 and causality explicitly, avoid per-frame find_f0 inside iaif
options.f0 = global_f0;  % used unless overwritten per frame (adaptive mode)
options.causality = 0;   % 0 -> 'noncausal' integration (backward), 1 -> 'causal'
options.remove_real_poles = 1; % Force removal of spectral tilt from VT model
% options.rho = 0.99;    % (optional) integrator leak, default is 0.99
% options.winfunc = @hamming;  % already default in iaif.m
%%% END NEW

% 6. Process each frame
fprintf('Processing frames with IAIF...\n');
results = cell(num_frames, 1);

for i = 1:num_frames
    if mod(i, 10) == 0
        fprintf('  Frame %d/%d\n', i, num_frames);
    end
    
    frame_data = frames{i};
    
    %%% NEW: if you use adaptive mode, update f0 per frame
    if strcmp(analysis_mode, 'adaptive') && isfield(pitch_info, 'f0') ...
            && ~isempty(pitch_info.f0)
        % You must define what pitch_info.f0(i) means (per frame)
        options.f0 = pitch_info.f0(i);
    end
    %%% END NEW
    
    % Apply IAIF to frame
    try
        [g, Hvt, e_ar, Hg] = iaif(frame_data, fs, options);
        
        if isa(g, 'signal')
            g = g.s;
        end
        
        % Store results
        results{i}.g = g(:);        % Glottal flow (IAIF)
        results{i}.Hvt = Hvt;       % Vocal tract filter
        results{i}.e_ar = e_ar(:);  % AR residual / error
        results{i}.Hg = Hg;         % Glottal filter
        results{i}.frame_idx = frame_indices(i);
        results{i}.frame_length = length(frame_data);
        results{i}.success = true;
    catch ME
        warning('Frame %d failed: %s', i, ME.message);
        results{i}.success = false;
        results{i}.g = [];
    end
end

fprintf('Done processing!\n');

% 8. Compare with ground truth
% NOTE: For a real comparison with glottal flow, this path should usually
% point to the *ug* file, e.g. '...-ug-....wav', not the speech signal.
gt_file = 'data/Glottal_signals_db/aa-ug-105Hz-8kHz.wav';  % Adjust path
compare_with_ground_truth(x, fs, results, frame_indices, 'IAIF', gt_file);