function [outputArg] = FIM_interval_complete_ob(par, data)
% Observed Fisher information matrix using complete likelihood function.
% usage: [outputArg] = FIM_interval_complete_ob(par, data)
% 
% arguments: (input)
%  par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'starttime', 'endtime', 'count', 'censored', 
%         'failure' and 'defective'.
% 
% arguments: (output)
%  outputArg - Observed Fisher information matrix using complete likelihood function.
%
% Example usage:
%  data = readtable('...\Data\CB_data.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:9) = 1; data.failure(13:17) = 2;
%  data.defective = zeros(height(data), 1);
%  par = [0.01, 138, 0.3, 4683, 3];
%  Result_complete_ob = FIM_interval_complete_ob(par, data)
%
    func_c = @(x) Obj_Ic(x, data, par);
    outputArg = hessian(func_c, par);
end

function [final_value] = Obj_Ic(initial_par, data, par_r)
    %% Parameters
    % theta
    p1 = initial_par(1); 
    a1 = initial_par(2); b1 = initial_par(3); 
    a2 = initial_par(4); b2 = initial_par(5);
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
            if ((fail_idx==1)&&(defective_idx==1)) % 
                L1 = count*log( value_g1 );
            elseif (fail_idx==2)&&(defective_idx==1)
                L1 = count*log( value_g2 );
            elseif (fail_idx==2)&&(defective_idx==2)
                L1 = count*log( value_g3 );
            end
        elseif (c1==1)&&(defective_idx>0)
            if (defective_idx==1)
                L1 = count*log( value_h1 );
            else
                L1= count*log( value_h2 );
            end
        elseif (c1==0)&&((fail_idx==0)||(defective_idx==0))
            if fail_idx==1
                L1 = count*log( value_g1 );
            elseif fail_idx==2
                prob_0 = value_g2_B/(value_g2_B+value_g3_B);
                L1 = count*( prob_0*log(value_g2) + (1-prob_0)*log(value_g3) );
            else
                prob_1 = value_g1_B/(value_g1_B+value_g2_B+value_g3_B);
                prob_2 = value_g2_B/(value_g1_B+value_g2_B+value_g3_B);
                L1 = count*( prob_1*log(value_g1) + prob_2*log(value_g2) + (1-prob_1-prob_2)*log(value_g3) );
            end
        elseif (c1==1)&&(defective_idx==0)
            prob_c = value_h1_B/(value_h1_B+value_h2_B);
            L1 = count*( prob_c*log(value_h1) + (1-prob_c)*log(value_h2) );
        end
        log_Qfunc = log_Qfunc + L1;
    end
    final_value = -log_Qfunc;

end

