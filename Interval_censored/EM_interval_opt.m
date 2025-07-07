function [par_em, iter_result] = EM_interval_opt(data, initial_par)
% EM algorithm: parameter estimation using EM algorithm with complete likelihood function.
% usage: [par_em, iter_result] = EM_interval_opt(data, initial_par)
% 
% arguments: (input)
%  data - Table with columns 'endtime', 'censored', 'failure' and 'defective'.
%  initial_par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%
% arguments: (output)
%  par_em - parameter estimation result.
%
%  iter_result - history of EM algorithm.
%
% Example usage:
%  data = readtable('...\Data\CB_data.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:9) = 1; data.failure(13:17) = 2;
%  data.defective = zeros(height(data), 1);
%  rng(7777777)
%  initial_par = rand(5, 1);
%  [par_em, iter_result] = EM_interval_opt(data, initial_par)
%
%
    iter = 100; i=1;
    par_iter = initial_par;
    iter_result = zeros(iter, 5);
    while i<=iter
        % Optimize
        par_iter(1) = value_pi;
        A = []; b = []; % A*x<=b
        Aeq = []; beq = []; % Aeq*x=beq 
        lb = [0, 0, 0, 0, 0]; 
        ub = [1, Inf, Inf, Inf, Inf]; % lb<=x<=ub
        options = optimset('MaxIter', 5000000, 'MaxFunEvals', 100000,'TolX',1e-10,'TolFun',1e-10);
        par_iter = fmincon(@(x) EM_interval_opt_func(x, data, par_iter), par_iter, A, b, Aeq, beq, lb, ub, [], options);
        iter_result(i, :) = par_iter;
        % Check stopping criteria
        if i>=6
            shift_sum = sum(abs(diff(iter_result(i-5:i, :), 1)), "all");
            wei_sum = all(sum(abs(diff(iter_result(i-5:i, :), 1)), 1)<(iter_result(i, :)./1000));
            if (shift_sum<0.005)||(wei_sum==1)     
                par_em = iter_result(i, :);
                break
            elseif i==iter
                par_em = [-1 -1 -1 -1 -1];
            end
        end
        i = i+1;
    end
        

end















