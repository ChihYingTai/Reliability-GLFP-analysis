function [par_em, iter_result] = EM_opt(data, initial_par)
% EM algorithm: parameter estimation using EM algorithm with complete likelihood function.
% usage: [par_em, iter_result] = EM_formula_opt(data, initial_par)
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
%  data = readtable('...\Data\Gateoxide.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:12) = 1; data.failure(20:44) = 2;
%  data.defective = zeros(height(data), 1);
%  rng(7777777)
%  initial_par = rand(5, 1);
%  [par_em, iter_result] = EM_opt(data, initial_par)
%
    iter = 100; i=1;
    par_iter = initial_par;
    iter_result = zeros(iter, 5);
    while i<=iter
        % Optimize
        value_pi = opt_pi(data, par_iter);
        par_iter(1) = value_pi;
        A = []; b = []; % A*x<=b
        Aeq = []; beq = []; % Aeq*x=beq 
        lb = [value_pi, 0, 0, 0, 0]; 
        ub = [value_pi, Inf, Inf, Inf, Inf]; % lb<=x<=ub
        options = optimset('MaxIter', 5000000, 'MaxFunEvals', 100000,'TolX',1e-10,'TolFun',1e-10);
        par_iter = fmincon(@(x) EM_opt_func(x, data, par_iter), par_iter, A, b, Aeq, beq, lb, ub, [], options);
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


function [value_pi] = opt_pi(data, par_r)
    B_p1 = par_r(1); 
    B_a1 = par_r(2); B_b1 = par_r(3); 
    B_a2 = par_r(4); B_b2 = par_r(5);
    %
    func_g1_B = @(x) B_p1.*wblpdf(x, B_a1, B_b1).*(1-wblcdf(x, B_a2, B_b2));
    func_g2_B = @(x) B_p1.*wblpdf(x, B_a2, B_b2).*(1-wblcdf(x, B_a1, B_b1));
    func_g3_B = @(x) (1-B_p1).*wblpdf(x, B_a2, B_b2);
    func_h1_B = @(x) B_p1.*(1-wblcdf(x, B_a1, B_b1)).*(1-wblcdf(x, B_a2, B_b2));
    func_h2_B = @(x) (1-B_p1).*(1-wblcdf(x, B_a2, B_b2));
    %
    fun_pi1 = @(x) func_g2_B(x) ./ ( func_g2_B(x)+func_g3_B(x) );
    fun_pi2 = @(x) func_g1_B(x) ./ ( func_g1_B(x)+func_g2_B(x)+func_g3_B(x) );
    fun_pi3 = @(x) func_h1_B(x) ./ ( func_h1_B(x)+func_h2_B(x) );
    fun_M = @(x) fun_pi1(x) + fun_pi2(x) - fun_pi1(x).*fun_pi2(x);
    %
    data_time = data.endtime; 
    data_a1_noc = data_time(((data.censored==0)&(data.failure==1)));
    data_a2_noc = data_time(((data.censored==0)&(data.failure==2)));
    data_m1_noc = data_time(((data.censored==0)&(data.failure==0)));
    data_cen = data_time((data.censored==1));
    %
    value_pi = (length(data_a1_noc) + sum(fun_pi3(data_cen)) + sum(fun_pi2(data_a2_noc)) + sum(fun_M(data_m1_noc)))/ length(data_time);

end















