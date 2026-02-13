classdef adc_model < handle
% ADC_MODEL - A class to model the behavior of a TI SAR ADC
% 28GSps 32-ways time interleaved ADC
% 7-bit resolution + 1-bit redundancy
% cdac weights: [64 32 16 8 8 4 2 1]
% calibration: gain, os, timing skew
    properties
        % parameter 
        fs = 28e9; % sampling frequncy
        N_ch = 32; % number of channels
        N_phase1 = 4; % number of phases in first stage
        Vref = 0.6; % reference votage

        % mismatch 
        % cof_gain = 0.01;
        % cof_os = 0.02;
        % cof_skew = 2e-12;
        % cof_cap = 0.016;

        cof_gain = 0;
        cof_os = 0;
        cof_skew = 0;
        cof_cap = 0.016;
        % core parameters
        cap_weights_single = [32, 16, 8, 4, 4, 2, 1]; % phsical capacitor weight
        dig_weights_diff = [64, 32, 16, 8, 8, 4, 2, 1]; % digital output weight (differential)
        
        % center code
        code_center; % ideal center code for differential output
        % capacitor mismatch
        cap_arr_pt; % [32*7] positive capacitor array including mismatch of top
        cap_arr_pb; % [32*7] positive capacitor array including mismatch of bottom
        cap_arr_nt; % [32*7] negative capacitor array including mismatch of top
        cap_arr_nb; % [32*7] negative capacitor array including mismatch of bottom

        % os
        os_comparator; % [32*1] offset of comparators
        % gain
        gain_ch; % [32*1] channel gain including afe only
        % skew
        skew_ch; % [4*1] channel skew

        % calibration registers
        reg_os_cal; % [32*1] offset calibration register
        reg_gain_cal; % [32*1] gain calibration register control calibration capacitor
        reg_dcdl;  % [4*1] delay calibration register

        % circuit parameters
        cunit = 5e-16; % unit capacitor size 0.5fF
        cp_top; % top plate parasitic capacitance
        cunit_cal = 2e-16; % unit capacitor size for calibration 0.1fF
        dcdl_step = 50e-15; % delay step of each dcdl tap 50fs
    end

    methods
        % adc initialization
        function obj = model_init(obj)
            fprintf("========Initializing ADC model========\n");
            obj.code_center = ceil(sum(obj.dig_weights_diff) / 2); % ideal center code for differential output
            obj.cp_top = 1/2*sum(obj.cap_weights_single)*obj.cunit; % top plate parasitic capacitance
            obj.init_hardware();
        end
        
        function init_hardware(obj)
            obj.cap_arr_pt = zeros(obj.N_ch, 7);
            obj.cap_arr_pb = zeros(obj.N_ch, 7);
            obj.cap_arr_nt = zeros(obj.N_ch, 7);
            obj.cap_arr_nb = zeros(obj.N_ch, 7);
            
            ideal_caps = obj.cap_weights_single * obj.cunit;  % ideal capacitor values
            % mismatch generation
            for ch = 1:obj.N_ch
                obj.cap_arr_nt(ch, :) = ideal_caps .* (1+obj.cof_cap*randn(1,7)); % 1.6% mismatch
                obj.cap_arr_nb(ch, :) = ideal_caps .* (1+obj.cof_cap*randn(1,7)); % 1.6% mismatch
                obj.cap_arr_pt(ch, :) = ideal_caps .* (1+obj.cof_cap*randn(1,7)); % 1.6% mismatch
                obj.cap_arr_pb(ch, :) = ideal_caps .* (1+obj.cof_cap*randn(1,7)); % 1.6% mismatch
            end 
            
            % os
            obj.os_comparator = obj.cof_os * randn(1,obj.N_ch); % 20mV std comparator offset
            % gain
            obj.gain_ch = 1 + obj.cof_gain * randn(1, obj.N_ch); % 1% std channel gain   
            % skew
            raw_skew = obj.cof_skew * randn(1, obj.N_phase1); % 2ps std channel skew
            obj.skew_ch = raw_skew - mean(raw_skew); % remove mean skew
            % registers reset
            obj.reg_os_cal = zeros(1, obj.N_ch);
            obj.reg_gain_cal = zeros(1, obj.N_ch);
            obj.reg_dcdl = zeros(1, obj.N_phase1);
        end

        function [dout_raw, t_act] = digitize_batch(obj, t_ideal, vin_func)  % vin_func: function handle generating diff_vin at given time
            N = length(t_ideal);  % sampling points
            dout_raw = zeros(1, N);
            t_act = zeros(1,N); % actual sampling time after skew and dcdl
            % get channel and phase index for each sample
            ch_list = mod(0:N-1, obj.N_ch) + 1; % channel index for each sample   
            ph_list = mod(ch_list - 1, obj.N_phase1) + 1; % phase index for each sample

            vcm = 0.5; % common mode voltage
            for i = 1:N
                ch = ch_list(i);
                ph = ph_list(i);
                % sampling
                current_skew = obj.skew_ch(ph) + obj.reg_dcdl(ph) * obj.dcdl_step;
                t_sample = t_ideal(i) + current_skew;
                t_act(i) = t_sample;
                % get vin
                vin_diff = vin_func(t_sample) * obj.gain_ch(ch);   % include afe gain

                % get capacitor array for this channel
                cap_pt = obj.cap_arr_pt(ch, :);
                cap_pb = obj.cap_arr_pb(ch, :);
                cap_nt = obj.cap_arr_nt(ch, :);
                cap_nb = obj.cap_arr_nb(ch, :);

                % capacitor for calibration
                cap_cal = obj.reg_gain_cal(ch) * obj.cunit_cal;  % for both sides

                % calculator the total capacitance
                c_tot_p = sum(cap_pt) + sum(cap_pb) + 2 * obj.cp_top + cap_cal;
                c_tot_n = sum(cap_nt) + sum(cap_nb) + 2 * obj.cp_top + cap_cal;

                % calculator the voltage sampled on the capacitor array
                v_sample_p = vin_diff/2 + vcm + obj.os_comparator(ch) / 2;
                v_sample_n = -vin_diff/2 +vcm - obj.os_comparator(ch) / 2;

                decisions  = zeros(1,8); % store comparator decisions

                v_node_p = v_sample_p;
                v_node_n = v_sample_n;
                % SAR conversion of cycles 1 to 7
                for bit = 1:7
                    if (v_node_p - v_node_n) >= 0
                        decisions(bit) = 1;
                        v_node_p = v_node_p - obj.Vref * cap_pt(bit) / c_tot_p;
                        v_node_n = v_node_n + obj.Vref * cap_nt(bit) / c_tot_n;
                    else 
                        decisions(bit) = 0;
                        v_node_p = v_node_p + obj.Vref * cap_pb(bit) / c_tot_p;
                        v_node_n = v_node_n - obj.Vref * cap_nb(bit) / c_tot_n;
                    end
                end
                % SAR conversion of cycle 8
                decisions(8) = (v_node_p - v_node_n) >= 0;

                % calculate digital output
                dout_raw(i) = sum(decisions .* obj.dig_weights_diff); 
            end 
        end
        % Calibration with batches
        function calib_with_batches(obj, N_total, batch_size, type, varargin)
            % prepare batches
            if strcmp(type, 'timing')
                fin = 12.99e9;
                amp = 0.35;
            else
                fin = 181.5e6;
                amp = 0.3;
            end

            vin_func = @(t) amp * sin(2*pi*fin*t);
            ts = 1/obj.fs;

            fprintf("========%s calibration with batches========\n",type);

            % OS calibration
            if strcmp(type, 'offset')
                acc_sum = zeros(1, obj.N_ch);
                acc_cnt = zeros(1, obj.N_ch);
                N_batches = ceil(N_total /batch_size);

                for b = 1:N_batches
                    t = ((b-1) * batch_size : b*batch_size - 1) * ts;
                    [d_raw, ~] = obj.digitize_batch(t, vin_func);

                    % calibration for current batch
                    for ch = 1:obj.N_ch
                        raw_ch = d_raw(ch: obj.N_ch :end);   % data for this channel
                        current_val = raw_ch - obj.code_center - obj.reg_os_cal(ch);
                        
                        acc_sum(ch) = acc_sum(ch) + sum(current_val);
                        acc_cnt(ch) = acc_cnt(ch) + length(current_val);
                    end
                end
                % update os registers
                average_os = acc_sum ./ acc_cnt;
                %obj.reg_os_cal = obj.reg_os_cal + 0.8 * average_os; % step size 0.4
                obj.reg_os_cal = average_os;
            
            elseif strcmp(type, 'gain')
                N_iter = varargin{1};
                for iter = 1 : N_iter
                    mae_sum = zeros(1,obj.N_ch);
                    mae_cnt = zeros(1,obj.N_ch);
                    N_batches = ceil(N_total / batch_size);

                    for b = 1: N_batches
                        t = ((b-1)*batch_size : b*batch_size - 1) * ts;
                        [d_raw, ~] = obj.digitize_batch(t, vin_func);

                        for ch = 1:obj.N_ch
                            raw_ch = d_raw(ch: obj.N_ch: end);
                            d_corr = raw_ch -obj.reg_os_cal(ch)-obj.code_center;

                            mae_sum(ch) = mae_sum(ch) +sum(abs(d_corr));   
                            mae_cnt(ch) = mae_cnt(ch) + length(d_corr);
                        end
                    end
                    % update gain registers
                    mae = mae_sum ./mae_cnt;
                    [ref_mae, ~] = min(mae);  % calculate the min mae

                    tol = 0.002 * ref_mae;
                    for ch = 1:obj.N_ch
                        if mae(ch) > (ref_mae + tol)
                            obj.reg_gain_cal(ch) = min([31, obj.reg_gain_cal(ch) + 1]); % max gain is 15
                        elseif mae(ch) < (ref_mae - tol)
                            obj.reg_gain_cal(ch) = max([0, obj.reg_gain_cal(ch) - 1]);
                        end   
                    end
                end
            
            % timing calibration
            elseif strcmp(type, 'timing')
                all_data = [];
                N_batches = ceil(N_total / batch_size);

                for b = 1:N_batches
                    t = ((b-1)*batch_size : b*batch_size - 1) * ts;
                    [d_raw, ~] = obj.digitize_batch(t, vin_func);
                    
                    d_corr = zeros(1, length(d_raw));
                    for ch = 1:obj.N_ch
                        idx = ch:obj.N_ch:length(d_raw);
                        d_corr(idx) = d_raw(idx) - obj.reg_os_cal(ch) - obj.code_center;
                    end 
                    all_data = [all_data, d_corr];
                end

                % calculate timing error for each phase
                ph_data = cell(1,obj.N_phase1);
                for p = 1:obj.N_phase1
                    ph_data{p} = all_data(p:obj.N_phase1:end); 
                end
                % calculate timing error
                for p = 1:obj.N_phase1
                    p_prev = mod(p-2, 4) + 1; p_next = mod(p, 4) + 1; % calculate the previous and next phase index
                    vec_c = ph_data{p};
                    vec_p = ph_data{p_prev};
                    vec_n = ph_data{p_next};
                    L = min([length(vec_c), length(vec_p), length(vec_n)])-1;
                    % extract equal length data
                    v_c = vec_c(1:L);
                    v_p = vec_p(1:L);
                    v_n = vec_n(1:L);

                    r_back = mean(v_p .* v_c);  % correlation with previous phase
                    r_fwd = mean(v_n .* v_c);   % correlation with next phase

                    diff = r_fwd - r_back;
                    threshold = 0.01 * abs(r_back);
                    if diff > threshold
                        obj.reg_dcdl(p) = max([0, obj.reg_dcdl(p) - 1]); % decrease delay
                    elseif diff < -threshold
                        obj.reg_dcdl(p) = min([15, obj.reg_dcdl(p) + 1]); % increase delay
                    end
                end
                fprintf('Timing Codes: [%d %d %d %d]\n', obj.reg_dcdl(1), obj.reg_dcdl(2), obj.reg_dcdl(3), obj.reg_dcdl(4));
            end
        end

        function run_analysis(obj)
            fprintf("========Running ADC analysis========\n");

            % parameter for analysis
            N_cal = 2^24;
            Batch_size = 2^16;

            N_test = 2^16;

            % reset and test
            fprintf("========Initial performance========\n");
            obj.reg_os_cal = zeros(1, obj.N_ch);
            obj.reg_gain_cal = zeros(1, obj.N_ch);
            obj.reg_dcdl = zeros(1, obj.N_phase1);
            fprintf("========Testing without calibration========\n");
            res0 = obj.test_performance(N_test, 'No calibration');

            % % os calibration
            % obj.calib_with_batches(N_cal, Batch_size, 'offset');
            % fprintf("========Testing after OS calibration========\n");
            % res_os = obj.test_performance(N_test, 'With OS calibration');

            % % gain calibration
            % obj.calib_with_batches(N_cal, Batch_size, 'gain', 15);
            % fprintf("========Testing after Gain calibration========\n");
            % res_gain = obj.test_performance(N_test, 'With Gain calibration');
            % 
            % % timing calibration
            % obj.calib_with_batches(N_cal, Batch_size, 'timing');
            % fprintf("========Testing after Timing calibration========\n");
            % res_timing = obj.test_performance(N_test, 'With Timing calibration');
            % 
            % % plot
            obj.plot_results({res0});
        end
    end

    methods (Access = private)
        function res = test_performance(obj, N, label)
            fin = 1.123e9; 
            amp = 0.39;
            vin_func = @(t) amp*sin(2*pi*fin*t);
            t = (0:N-1) / obj.fs;
            % adc digitization
            [d_raw, ~] = obj.digitize_batch(t, vin_func);
            % data correction digital domain
            d_corr = zeros(1,length(d_raw));
            for ch = 1:obj.N_ch
                d_corr(ch:obj.N_ch:end) = d_raw(ch:obj.N_ch:end) - obj.reg_os_cal(ch) - obj.code_center;
            end
            % fft
            w = blackman(N)';
            d = (d_corr - mean(d_corr)) .* w;
            D = fft(d);
            D = D(1:N/2);
            P = abs(D).^2;
            P = P / sum(w.^2);
            [~, i] = max(P(2:end));
            i = i + 1;
            if i > 5
                % signal bins
                sig_bins = (i-3):(i+3);
                sig_bins = sig_bins(sig_bins>1 & sig_bins<=N/2);
                Psig = sum(P(sig_bins));
                 % noise bins (exclude signal bins + DC)
                noise_bins = 2:N/2;
                noise_bins = setdiff(noise_bins, sig_bins);
                Pnoise = sum(P(noise_bins));
                SNR = 10*log10(Psig/Pnoise);
                enob = (SNR - 1.76)/6.02;
            else
                enob = 0;
            end
            spec = 10*log10(P + eps);
            spec = spec - max(spec);

            res.label = label;
            res.spec = spec;
            res.enob = enob;
            res.f = linspace(0, obj.fs/2,N/2);

            % plot waveform
            N_plot = min(512,N);
            scale = obj.Vref *2 / (3*obj.code_center);
            res.t_ns = t(1:N_plot) * 1e9;  % ns time scale
            res.v_in = vin_func(t(1:N_plot));
            res.v_adc = d_corr(1:N_plot) * scale;  % convert code to voltage
        end

        function res = test_performance2(obj, N, label)
            % 1. [关键修改] 计算相干采样频率
            % 目标频率附近寻找
            target_fin = 1.123e9; 
            
            % 计算周期数 M (四舍五入到最近的整数)
            M = round(target_fin / obj.fs * N);
            
            % 强制 M 为奇数 (确保与 N=2^k 互质, 防止重复相位)
            if mod(M, 2) == 0
                M = M + 1;
            end
            
            % 反推精确的输入频率
            fin = (M / N) * obj.fs;
            
            % 2. 生成信号
            amp = 0.395; % 0.9Vpp
            vin_func = @(t) amp * sin(2*pi*fin*t);
            t = (0:N-1) / obj.fs;
            
            % 3. ADC 采样
            [d_raw, ~] = obj.digitize_batch(t, vin_func);
            
            % 4. 数字处理
            d_corr = zeros(1, length(d_raw));
            for ch = 1:obj.N_ch
                % 注意：这里修正了变量名和减法逻辑
                d_corr(ch:obj.N_ch:end) = d_raw(ch:obj.N_ch:end) - obj.reg_os_cal(ch) - obj.code_center;
            end
            
            % 5. [关键修改] FFT 分析 (相干采样不需要加窗，或者用矩形窗)
            % win = blackman(N)'; % <-- 删掉或注释掉 Blackman 窗
            win = ones(1, N);     % <-- 改用矩形窗 (无窗)
            
            data_ac = d_corr - mean(d_corr);
            
            % 计算频谱
            spec_mag = abs(fft(data_ac .* win));
            spec = 20*log10(spec_mag(1:N/2) + eps);
            spec = spec - max(spec); % 归一化到 0dBFS
            
            % 6. 计算指标 (SNDR/ENOB)
            [~, idx] = max(spec);
            
            % 信号功率 (因为是相干采样，能量高度集中，通常只取峰值点即可，或者取临近3点以防万一)
            % 但为了稳健，还是保留邻域求和
            if idx > 5
                sig_p = sum(10.^(spec(idx-3:idx+3)/10)); 
                noi_p = sum(10.^(spec(5:end)/10)) - sig_p;
                
                if noi_p > 0
                    sndr = 10*log10(sig_p/noi_p);
                    enob = (sndr - 1.76) / 6.02;
                else
                    enob = 0; % 避免 log(0)
                end
            else
                enob = 0;
            end

            res.label = label;
            res.spec = spec;
            res.enob = enob;
            res.f = linspace(0, obj.fs/2, N/2);

            % 绘图波形数据
            N_plot = min(500, N);
            scale = obj.Vref *2 / (3*obj.code_center);
            res.t_ns = t(1:N_plot) * 1e9;  
            res.v_in = vin_func(t(1:N_plot));
            res.v_adc = d_corr(1:N_plot) * scale; 
        end
        function plot_results(obj, res_list)
            figure('Name', 'Analysis summary', 'Position', [50, 50, 1200, 500], 'Color', 'black');
            for i = 1:length(res_list)
                r = res_list{i};
                % spectrum
                subplot(length(res_list),2,2*i-1);
                plot(r.f/1e9, r.spec, 'b', 'LineWidth', 0.8); grid on;
                ylim([-100 0]);
                xlim([0 obj.fs/2e9]);
                title(sprintf('%s ENoB = %.2f bits',r.label, r.enob));
                if i == length(res_list)
                    xlabel('Frequency (GHz)');
                end
                ylabel('Magnitude (dB)');
                %waveform
                subplot(length(res_list),2,2*i);
                plot(r.t_ns, r.v_in, 'g', 'LineWidth',2); hold on;
                stairs(r.t_ns, r.v_adc, 'r--'); grid on;
                ylim([-0.6 0.6]);
                xlim([0 r.t_ns(end)/3]); % Zoom in
                title(sprintf('%s Input and ADC output waveforms', r.label));
                if i == length(res_list)
                    xlabel('Time (ns)');
                end
            end
        end
    end
end