function [outputArg] = FIM_interval_missing_ob(par, data)
% Observed Fisher information matrix using missing data.
% usage: [outputArg] = FIM_interval_missing_ob(par, data)
% 
% arguments: (input)
%  par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'starttime', 'endtime', 'count', 'censored', 
%         'failure' and 'defective'.
% 
% arguments: (output)
%  outputArg - Observed Fisher information matrix using missing data.
%
% Example usage:
%  data = readtable('...\Data\CB_data.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:9) = 1; data.failure(13:17) = 2;
%  data.defective = zeros(height(data), 1);
%  par = [0.01, 138, 0.3, 4683, 3];
%  Result_missing_ob = FIM_missing_ob(par, data)
%
    func_m = @(x) Obj_Im(x, data, par);
    outputArg = hessian(func_m, par);
end

function [final_value] = Obj_Im(par_r1, data, par_r)
    %% Parameters
    % theta
    p1 = par_r1(1); 
    a1 = par_r1(2); b1 = par_r1(3); 
    a2 = par_r1(4); b2 = par_r1(5);
    % theta^{(r)}
    B_p1 = par_r(1); 
    B_a1 = par_r(2); B_b1 = par_r(3); 
    B_a2 = par_r(4); B_b2 = par_r(5);
    %% Functions
    % Functions for theta^{(r)}
    func_g1_B = @(x) B_p1.*wblpdf(x, B_a1, B_b1).*(1-wblcdf(x, B_a2, B_b2));
    func_g2_B = @(x) B_p1.*wblpdf(x, B_a2, B_b2).*(1-wblcdf(x, B_a1, B_b1));
    % Functions for theta
    func_g1 = @(x) p1.*wblpdf(x, a1, b1).*(1-wblcdf(x, a2, b2));
    func_g2 = @(x) p1.*wblpdf(x, a2, b2).*(1-wblcdf(x, a1, b1));
    %%
    log_Qfunc = 0;
    for i = 1:height(data)
        d1 = [data.starttime(i) data.endtime(i)];
        count = data.count(i);
        c1 = data.censored(i); 
        defective_idx = data.defective(i); 
        fail_idx = data.failure(i);
        %
        value_g1 = integral(func_g1, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g2 = integral(func_g2, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g3 = (1-p1).*(wblcdf(d1(2), a2, b2)-wblcdf(d1(1), a2, b2));
        value_h1 = p1.*(1-wblcdf(d1(1), a1, b1)).*(1-wblcdf(d1(1), a2, b2));
        value_h2 = (1-p1).*(1-wblcdf(d1(1), a2, b2));
        %
        value_g1_B = integral(func_g1_B, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g2_B = integral(func_g2_B, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g3_B = (1-B_p1).*(wblcdf(d1(2), B_a2, B_b2)-wblcdf(d1(1), B_a2, B_b2));
        value_h1_B = B_p1.*(1-wblcdf(d1(1), B_a1, B_b1)).*(1-wblcdf(d1(1), B_a2, B_b2));
        value_h2_B = (1-B_p1).*(1-wblcdf(d1(1), B_a2, B_b2));
        %
        if (c1==0)&&(fail_idx>0)&&(defective_idx>0)
            continue
        elseif (c1==1)&&(defective_idx>0)
            continue
        elseif (c1==0)&&((fail_idx==0)||(defective_idx==0))
            value_1 = value_g1/(value_g1+value_g2+value_g3);
            value_2 = value_g2/(value_g1+value_g2+value_g3);
            value_c = value_h1/(value_h1+value_h2);
            if fail_idx==1
                continue
            elseif fail_idx==2
                prob_0_B = value_g2_B/(value_g2_B+value_g3_B);
                L1 = count*( prob_0_B*log(value_2) + (1-prob_0_B)*log(1-value_1-value_2) );
            else
                prob_1_B = value_g1_B/(value_g1_B+value_g2_B+value_g3_B);
                prob_2_B = value_g2_B/(value_g1_B+value_g2_B+value_g3_B);
                L1 = count*( prob_1_B*log(value_1) + prob_2_B*log(value_2) + (1-prob_1_B-prob_2_B)*log(1-value_1-value_2) );
            end
        elseif (c1==1)&&(defective_idx==0)
            prob_c = value_h1_B/(value_h1_B+value_h2_B);
            L1 = count*( prob_c*log(value_c) + (1-prob_c)*log(1-value_c) );
        end
        log_Qfunc = log_Qfunc + L1;
    end
    final_value = log_Qfunc;

end

















