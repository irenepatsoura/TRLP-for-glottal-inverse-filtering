% Script to compute NAQ, QOQ, HRF, and H1-H2 metrics
% Iterates over vowels /a/, /e/, /i/ and f0 groups.

addpath("core","eval","framework","pipeline")

vowels = {'aa', 'eh', 'ih'};
vowel_names = {'/a/', '/e/', '/i/'};

f0_groups = {
    'Male', [105, 115, 130, 145];
    'Female', [205, 210, 230, 255, 310, 410]
};

methods = {'QCP', 'IAIF', 'TRLP'}; % Methods to compare against GT
metrics_list = {'NAQ', 'QOQ', 'HRF', 'H1H2'};

% Common parameters
fs_target = 8000;
frame_len_ms = 50;
frame_shift_ms = 25;

fprintf('Starting Glottal Metrics Analysis...\n');

% Data structure to hold all results for plotting
% all_results.(metric).(group).(vowel).(method) = [errors]
all_results = struct();
for m = 1:length(metrics_list)
    all_results.(metrics_list{m}) = struct();
end

for v_idx = 1:length(vowels)
    v_code = vowels{v_idx};
    v_name = vowel_names{v_idx};
    
    for g_idx = 1:size(f0_groups, 1)
        group_name = f0_groups{g_idx, 1};
        f0_list = f0_groups{g_idx, 2};
        
        % Initialize accumulators for this vowel-group
        for m = 1:length(metrics_list)
            metric = metrics_list{m};
            for meth = 1:length(methods)
                all_results.(metric).(group_name).(v_code).(methods{meth}) = [];
            end
        end
        
        for f0 = f0_list
            speech_file = sprintf('data/Glottal_signals_db/%s-pout-%dHz-8kHz.wav', v_code, f0);
            gt_file = sprintf('data/Glottal_signals_db/%s-ug-%dHz-8kHz.wav', v_code, f0);
            
            if ~exist(speech_file, 'file') || ~exist(gt_file, 'file')
                continue;
            end
            
            fprintf('Processing %s %s (%d Hz)...\n', group_name, v_name, f0);
            
            % Load signals
            [x, fs] = audioread(speech_file);
            if size(x, 2) > 1, x = mean(x, 2); end
            x = x(:);
            
            [gt, fs_gt] = audioread(gt_file);
            if size(gt, 2) > 1, gt = mean(gt, 2); end
            gt = gt(:);
            
            if fs ~= fs_target, x = resample(x, fs_target, fs); fs = fs_target; end
            if fs_gt ~= fs_target, gt = resample(gt, fs_target, fs_gt); end
            
            if skewness(x) < 0, x = -x; end
            [b_hp, a_hp] = butter(2, 50/(fs/2), 'high');
            x = filtfilt(b_hp, a_hp, x);
            
            frame_len = round(frame_len_ms * fs / 1000);
            frame_shift = round(frame_shift_ms * fs / 1000);
            [frames, frame_indices] = create_fixed_frames(x, frame_len, frame_shift);
            [gt_frames, ~] = create_fixed_frames(gt, frame_len, frame_shift);
            
            % Run Methods
            % QCP Options
            qcp_opts = struct('dq', 0.4, 'pq', 0.05, 'nramp', round(fs/8000*7), ...
                              'causality', 0, 'f0', f0, 'remove_real_poles', 1);
            
            % IAIF Options
            iaif_opts = struct('p', round(fs/1000)+2, 'g', 4, 'rho', 0.99, ...
                               'causality', 0, 'f0', f0, 'remove_real_poles', 1);
            
            % TRLP Initialization
            p_trlp = round(fs/1000) + 2;
            a_prev = zeros(p_trlp, 1);
            rho_trlp = 0.99;
            
            n_frames = min(length(frames), length(gt_frames));
            for i = 1:n_frames
                frame = frames{i};
                gt_frame = gt_frames{i};
                
                % GT Metrics
                m_gt = compute_glottal_metrics(gt_frame, fs, f0);
                
                % QCP
                % Pass as signal object to match run_qcp.m
                [g_qcp, ~] = qcp(signal(frame, fs), qcp_opts);
                m_qcp = compute_glottal_metrics(g_qcp, fs, f0);
                
                % IAIF
                [g_iaif, ~] = iaif(frame, fs, iaif_opts);
                m_iaif = compute_glottal_metrics(g_iaif, fs, f0);
                
                % TRLP
                a_trlp = mytrlp(frame, p_trlp, a_prev, 0.5);
                a_prev = a_trlp;
                res_trlp = filter([1; -a_trlp], 1, frame);
                g_trlp = filter(1, [1 -rho_trlp], res_trlp);
                m_trlp = compute_glottal_metrics(g_trlp, fs, f0);
                
                % Store Errors (Est - GT)
                for m = 1:length(metrics_list)
                    metric = metrics_list{m};
                    val_gt = m_gt.(metric);
                    
                    if ~isnan(val_gt)
                        % QCP
                        val_est = m_qcp.(metric);
                        if ~isnan(val_est)
                            all_results.(metric).(group_name).(v_code).QCP(end+1) = val_est - val_gt;
                        end
                        
                        % IAIF
                        val_est = m_iaif.(metric);
                        if ~isnan(val_est)
                            all_results.(metric).(group_name).(v_code).IAIF(end+1) = val_est - val_gt;
                        end
                        
                        % TRLP
                        val_est = m_trlp.(metric);
                        if ~isnan(val_est)
                            all_results.(metric).(group_name).(v_code).TRLP(end+1) = val_est - val_gt;
                        end
                    end
                end
            end
        end
    end
end

% Plotting
for m = 1:length(metrics_list)
    metric = metrics_list{m};
    figure('Name', [metric ' Error Analysis'], 'NumberTitle', 'off', 'Color', 'w');
    
    % Subplot 1: Male
    subplot(1, 2, 1);
    plot_group_metrics(all_results.(metric).Male, vowels, vowel_names, methods, ['Male - ' metric ' Error']);
    
    % Subplot 2: Female
    subplot(1, 2, 2);
    plot_group_metrics(all_results.(metric).Female, vowels, vowel_names, methods, ['Female - ' metric ' Error']);
end

% Print Tables
fprintf('\n========================================\n');
fprintf('       Error Analysis Results       \n');
fprintf('========================================\n');

for m = 1:length(metrics_list)
    metric = metrics_list{m};
    fprintf('\nMetric: %s (Error = Est - GT)\n', metric);
    fprintf('----------------------------------------\n');
    
    for g_idx = 1:size(f0_groups, 1)
        group_name = f0_groups{g_idx, 1};
        fprintf('Group: %s\n', group_name);
        
        % Header
        fprintf('%-10s', 'Vowel');
        for meth = 1:length(methods)
            fprintf('%-25s', methods{meth});
        end
        fprintf('\n');
        
        for v_idx = 1:length(vowels)
            v_code = vowels{v_idx};
            v_name = vowel_names{v_idx};
            fprintf('%-10s', v_name);
            
            for meth = 1:length(methods)
                method = methods{meth};
                data = all_results.(metric).(group_name).(v_code).(method);
                if ~isempty(data)
                    mu = mean(data);
                    sigma = std(data);
                    fprintf('%-25s', sprintf('%.4f +/- %.4f', mu, sigma));
                else
                    fprintf('%-25s', 'N/A');
                end
            end
            fprintf('\n');
        end
        fprintf('\n');
    end
end

function plot_group_metrics(group_data, vowels, vowel_names, methods, title_str)
    n_vowels = length(vowels);
    n_methods = length(methods);
    
    means = zeros(n_vowels, n_methods);
    stds = zeros(n_vowels, n_methods);
    
    for v = 1:n_vowels
        v_code = vowels{v};
        for meth = 1:n_methods
            method = methods{meth};
            data = group_data.(v_code).(method);
            if ~isempty(data)
                means(v, meth) = mean(data);
                stds(v, meth) = std(data);
            else
                means(v, meth) = 0;
                stds(v, meth) = 0;
            end
        end
    end
    
    b = bar(means);
    hold on;
    
    % Error bars
    ngroups = size(means, 1);
    nbars = size(means, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    
    for i = 1:nbars
        % Calculate center of each bar
        if isprop(b(i), 'XEndPoints')
             x = b(i).XEndPoints;
        else
             x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        end
        errorbar(x, means(:,i), stds(:,i), 'k', 'linestyle', 'none');
    end
    
    set(gca, 'XTickLabel', vowel_names);
    title(title_str);
    ylabel('Error (Est - GT)');
    legend(methods, 'Location', 'best');
    grid on;
    hold off;
end
