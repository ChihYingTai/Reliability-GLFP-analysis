function [data_table] = Data_func(info)
% Generate data with specified masking level.
% usage: [data_table] = Data_func(info)
% 
% arguments: (input)
%  info - Structure with n_sample: samplesize;
%                       censored_t: right censored time;
%                       mask_candi: masking levels;
%                       par: ture parameters.
% 
% arguments: (output)
%  data_table - Data table.
%
% Example usage:
%  info.n_sample = 100; 
%  info.censored_t = 185.5; 
%  info.mask_candi = [0, 0.5, 0.8];
%  info.par = [0.1, 1.5, 0.45, 183, 7];  
%  data_table = Data_func(info)
%
    %% Generate raw data
    par = info.par;
    n_sample = info.n_sample;
    censored_t = info.censored_t;
    mask_candi = info.mask_candi;
    p = par(1); a_1 = par(2); b_1 = par(3); a_2 = par(4); b_2 = par(5); 
    % 
    flag = 1;
    while flag==1
        defect_idx = mnrnd(1, [p 1-p], n_sample);
        flag = any(sum(defect_idx, 1)<5);
    end
    %
    life_1 = wblrnd(a_1, b_1, n_sample, 1); life_2 = wblrnd(a_2, b_2, n_sample, 1);
    generate_data = zeros(n_sample, 5);
    for i = 1:n_sample
        defect_1 = defect_idx(i, :);
        if defect_1(2)==1
            generate_data(i, :) = [life_2(i), 2, 2, life_1(i), life_2(i)];
        else
            [c_time, min_idx] = min([life_1(i), life_2(i)]);
            generate_data(i, :) = [c_time, min_idx, 1, life_1(i), life_2(i)];
        end
    end
    data_table = array2table(generate_data, 'VariableNames',{'endtime', 'failure_0', 'defective', 'life_1', 'life_2'});
    data_table.endtime_raw = data_table.endtime;
    data_table.failure_raw = data_table.failure_0;
    data_table.defective_raw = data_table.defective;
    [~, s1] = sort(data_table.endtime);
    data_table = data_table(s1, :);
    data_table.defective(1:end) = 0;
    % Add censored data and add mask to it
    idx_cen = find(data_table.endtime>censored_t);
    data_table.censored = zeros(height(data_table), 1);
    if length(idx_cen)>=1
        data_table.censored(idx_cen) = 1;
        data_table.endtime(idx_cen) = censored_t;
        data_table.failure_0(idx_cen) = 0;
    end
    %% Add mask to failure data
    for i = 1:length(mask_candi)
        mask_ratio = mask_candi(i);
        if mask_ratio==0
            continue
        end
        name_f = "failure_"+string(mask_ratio*10);
        data_table.(name_f) = data_table.failure_0;
        %
        idx_f1 = find(data_table.(name_f)==1);
        mask_nf1 = floor(length(idx_f1)*mask_ratio);
        if length(idx_f1)>2
            if length(idx_f1)-mask_nf1>2
                s1 = randsample(idx_f1(2:end-1), mask_nf1);
                data_table.(name_f)(s1) = 0;
            else
                data_table.(name_f)(idx_f1(2:end-1)) = 0;
            end
        end
        %
        idx_f2 = find(data_table.(name_f)==2);
        mask_nf2 = floor(length(idx_f2)*mask_ratio);
        if length(idx_f2)>2
            if length(idx_f2)-mask_nf2>2
                s2 = randsample(idx_f2(2:end-1), mask_nf2);
                data_table.(name_f)(s2) = 0;
            else
                data_table.(name_f)(idx_f2(2:end-1)) = 0;
            end
        end
    end




end









