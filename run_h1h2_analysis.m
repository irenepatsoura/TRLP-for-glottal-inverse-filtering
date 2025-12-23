% Script to compute H1-H2 metrics for Ground Truth, QCP, IAIF, and TRLP
% Iterates over vowels /a/, /e/, /i/, /u/ (mapped to aa, eh, ih, uh/uw?)
% and f0 groups.

addpath("core","eval","framework","pipeline")

% Define vowels and their file codes
% Based on file list: aa, ae, eh, ih. 'u' is missing from the list I saw?
% Let's assume 'aa', 'eh', 'ih' are available. I didn't see 'u' or 'uw'.
% I saw 'aa', 'ae', 'eh', 'ih'.
% The user asked for /a, e, i, ou/.
% 'aa' -> /a/
% 'eh' -> /e/
% 'ih' -> /i/
% 'uw' or 'uh' -> /ou/ (u). I need to check if 'uw' exists.
% Let's check for 'uw' or 'uh' specifically.

vowels = {'aa', 'eh', 'ih'}; % Add 'uw' if found
vowel_names = {'/a/', '/e/', '/i/'};

% Define f0 groups
% Male: [105, 115, 130, 145]
% Female: [205, 210, 230, 255, 310, 410]
f0_groups = {
    'Male', [105, 115, 130, 145];
    'Female', [205, 210, 230, 255, 310, 410]
};

methods = {'GT', 'QCP', 'IAIF', 'TRLP'};
results_table = {};

% Common parameters
fs_target = 8000;
frame_len_ms = 50;
frame_shift_ms = 25; % 50% overlap

fprintf('Starting H1-H2 Analysis...\n');

for v_idx = 1:length(vowels)
    v_code = vowels{v_idx};
    v_name = vowel_names{v_idx};
    
    for g_idx = 1:size(f0_groups, 1)
        group_name = f0_groups{g_idx, 1};
        f0_list = f0_groups{g_idx, 2};
        
        % Accumulators for this vowel-group combination
        h1h2_acc = struct('GT', [], 'QCP', [], 'IAIF', [], 'TRLP', []);
        
        for f0 = f0_list
            % Construct filename pattern
            % e.g., aa-pout-105Hz-8kHz.wav
            % We need the SPEECH signal for processing, and GT for reference.
            % Speech: [v_code]-pout-[f0]Hz-8kHz.wav
            % GT:     [v_code]-ug-[f0]Hz-8kHz.wav
            
            speech_file = sprintf('data/Glottal_signals_db/%s-pout-%dHz-8kHz.wav', v_code, f0);
            gt_file = sprintf('data/Glottal_signals_db/%s-ug-%dHz-8kHz.wav', v_code, f0);
            
            if ~exist(speech_file, 'file') || ~exist(gt_file, 'file')
                fprintf('Skipping missing file: %s\n', speech_file);
                continue;
            end
            
            fprintf('Processing %s (%d Hz)...\n', v_name, f0);
            
            % Load signals
            [x, fs] = audioread(speech_file);
            if size(x, 2) > 1, x = mean(x, 2); end
            x = x(:);
            
            [gt, fs_gt] = audioread(gt_file);
            if size(gt, 2) > 1, gt = mean(gt, 2); end
            gt = gt(:);
            
            % Resample if needed (though filenames say 8kHz)
            if fs ~= fs_target, x = resample(x, fs_target, fs); fs = fs_target; end
            if fs_gt ~= fs_target, gt = resample(gt, fs_target, fs_gt); end
            
            % Pre-processing (Highpass, Polarity)
            if skewness(x) < 0, x = -x; end
            [b_hp, a_hp] = butter(2, 50/(fs/2), 'high');
            x = filtfilt(b_hp, a_hp, x);
            
            % Create Frames
            frame_len = round(frame_len_ms * fs / 1000);
            frame_shift = round(frame_shift_ms * fs / 1000);
            [frames, frame_indices] = create_fixed_frames(x, frame_len, frame_shift);
            
            % --- Run Methods ---
            
            % 1. GT Frames
            % Extract GT frames corresponding to speech frames
            % (Assuming alignment is roughly correct, or we just analyze GT frames independently)
            % Better to analyze GT frames extracted from the GT file using same indices
            [gt_frames, ~] = create_fixed_frames(gt, frame_len, frame_shift);
            
            % 2. QCP
            qcp_opts = struct('dq', 0.4, 'pq', 0.05, 'nramp', round(fs/8000*7), 'causality', 0, 'f0', f0, 'remove_real_poles', 1);
            
            % 3. IAIF
            iaif_opts = struct('p', (fs/1000 + 2), 'g', 4, 'f0', f0, 'causality', 0, 'remove_real_poles', 1);
            
            % 4. TRLP
            trlp_p = fs/1000 + 2;
            trlp_rho = 0.99;
            a_prev = zeros(trlp_p, 1);
            
            for i = 1:length(frames)
                frm_speech = frames{i};
                if i > length(gt_frames), break; end
                frm_gt = gt_frames{i};
                
                % --- GT H1-H2 ---
                val = compute_h1h2(frm_gt, fs, f0);
                h1h2_acc.GT = [h1h2_acc.GT; val];
                
                % --- QCP ---
                try
                    frame_sig = signal(frm_speech, fs);
                    [sg, ~, ~, ~] = qcp(frame_sig, qcp_opts);
                    val = compute_h1h2(sg.s, fs, f0);
                    h1h2_acc.QCP = [h1h2_acc.QCP; val];
                catch
                end
                
                % --- IAIF ---
                try
                    [g_iaif, ~, ~, ~] = iaif(frm_speech, fs, iaif_opts);
                    if isa(g_iaif, 'signal'), g_iaif = g_iaif.s; end
                    val = compute_h1h2(g_iaif, fs, f0);
                    h1h2_acc.IAIF = [h1h2_acc.IAIF; val];
                catch
                end
                
                % --- TRLP ---
                try
                    a = mytrlp(frm_speech, trlp_p, a_prev, 0.5);
                    a_prev = a;
                    res = filter([1; -a], 1, frm_speech);
                    g_trlp = filter(1, [1 -trlp_rho], res);
                    val = compute_h1h2(g_trlp, fs, f0);
                    h1h2_acc.TRLP = [h1h2_acc.TRLP; val];
                catch
                end
            end
        end
        
        % Calculate Mean/Std for this group
        row = {v_name, group_name};
        for m = 1:length(methods)
            vals = h1h2_acc.(methods{m});
            if isempty(vals)
                str = 'N/A';
            else
                str = sprintf('%.2f (%.2f)', mean(vals), std(vals));
            end
            row = [row, {str}];
        end
        results_table = [results_table; row];
    end
end

% Display Table
fprintf('\n\n=== H1-H2 Metrics (Mean (Std) dB) ===\n');
fprintf('%-10s %-10s %-15s %-15s %-15s %-15s\n', 'Vowel', 'Group', 'GT', 'QCP', 'IAIF', 'TRLP');
for i = 1:size(results_table, 1)
    fprintf('%-10s %-10s %-15s %-15s %-15s %-15s\n', results_table{i, :});
end
