function [par_mle, output_result] = MLE_interval_opt(data, initial_par)
% MLE: parameter estimation using MLE with incomplete likelihood function.
% usage: [par_mle, output_result] = MLE_interval_opt(data, initial_par)
% 
% arguments: (input)
%  data - Table with columns 'endtime', 'censored', 'failure' and 'defective'.
% 
%  initial_par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%
% arguments: (output)
%  par_mle - parameter estimation result.
%
%  output_result - history of different initials.
%
% Example usage:
%  data = readtable('...\Data\CB_data.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:9) = 1; data.failure(13:17) = 2;
%  data.defective = zeros(height(data), 1);
%  rng(7777777)
%  initial_par = rand(5, 1);
%  [par_mle, output_result] = MLE_interval_opt(data, initial_par)
%
%
    iter_mle = 10;
    final_result = zeros(iter_mle, 16);
    i=1; j=0;
    A = []; b = []; % A*x<=b
    Aeq = []; beq = []; % Aeq*x=beq 
    lb = [0, 0, 0, 0, 0]; 
    ub = [1, Inf, Inf, Inf, Inf]; % lb<=x<=ub
%     options = optimset('MaxIter', 5000000, 'MaxFunEvals', 100000,'TolX',1e-10,'TolFun',1e-10);
    while i<=iter_mle
        try
            if (i==1)&&(j==0)
                x0 = initial_par;
                j = 1;
            else
                x0 = initial_par + [rand, rand*5, rand*5, rand*10, rand*5]';
            end
            [par_opt, fval, ~, ~, ~, ~, hessian] = fmincon(@(x) MLE_interval_opt_func(x, data), x0, A, b, Aeq, beq, lb, ub, []); % , options
            hessian_inv = inv(hessian);
            final_result(i, :) = [-fval par_opt' diag(hessian_inv)' x0'];
            i = i+1;
        end
    end
    [~, s1] = max(final_result(:, 1));
    output_result = final_result(s1, :);
    par_mle = final_result(s1, 2:6); 
end

















