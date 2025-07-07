function [result_tp1, result_tp2, result_systp] = FIM_SE_CI_tp_func(Cov_par, par, candicate_q)
% Calculate SE and CI for tp1, tp2 and system tp.
% usage: [result_tp1, result_tp2, result_systp] = FIM_CI_tp_func(Cov_par, par, candicate_q)
% 
% arguments: (input)
%  Cov_par - Matrix of covariance matrix.
%
%  par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%
%  candicate_q - Scalar of candicate quantile.
% 
% arguments: (output)
%  result_tp1 - the estimates, SE and CI of tp1.
%
%  result_tp2 - the estimates, SE and CI of tp2.
%
%  result_systp - the estimates, SE and CI of system tp.
%
% Example usage:
%  data = readtable('...\Data\Gateoxide.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:12) = 1; data.failure(20:44) = 2;
%  data.defective = zeros(height(data), 1);
%  par = [0.4, 1.3, 0.1, 185, 8];
%  Cov_par = FIM_complete_ob(par, data);
%  [result_tp1, result_tp2, result_systp] = FIM_CI_tp_func(Cov_par, par, candicate_q)
%
    % Parameters
    p = par(1);
    a1 = par(2); b1 = par(3);
    a2 = par(4); b2 = par(5);
    % tp1
    tp1_hat = a1.*(-log(1-p)).^(1/b1);
    v1 = (-log(1-p)).^(1/b1);
    v2 = -a1.*(-log(1-p)).^(1/b1).*b1.^(-2).*log(-log(1-p));
    se_tp1 = sqrt( v1.^(2).*Cov_par(2,2) + v2.^(2).*Cov_par(3,3) +2.*v1.*v2.*Cov_par(2,3) );
    CI_U_tp1 = tp1_hat.*exp(1.96.*se_tp1./tp1_hat);
    CI_L_tp1 = tp1_hat.*exp(-1.96.*se_tp1./tp1_hat);
    CI_tp1 = zeros(1, 2*length(p));
    CI_tp1(1:2:(2*length(p)-1)) = CI_L_tp1;
    CI_tp1(2:2:(2*length(p))) = CI_U_tp1;
    % tp2
    tp2_hat = a2.*(-log(1-p)).^(1/b2);
    v1 = (-log(1-p)).^(1/b2);
    v2 = -a2.*(-log(1-p)).^(1/b2).*b2.^(-2).*log(-log(1-p));
    se_tp2 = sqrt( v1.^(2).*Cov_par(4,4) + v2.^(2).*Cov_par(5,5) +2.*v1.*v2.*Cov_par(4,5) );
    CI_U_tp2 = tp2_hat.*exp(1.96.*se_tp2./tp2_hat);
    CI_L_tp2 = tp2_hat.*exp(-1.96.*se_tp2./tp2_hat);
    CI_tp2 = zeros(1, 2*length(p));
    CI_tp2(1:2:(2*length(p)-1)) = CI_L_tp2;
    CI_tp2(2:2:(2*length(p))) = CI_U_tp2;
    % tp - system
    [systp_hat, se_systp] = system_tp_simulation(Cov_par, par, candicate_q);
    CI_U_systp = systp_hat.*exp(1.96.*se_systp./systp_hat);
    CI_L_systp = systp_hat.*exp(-1.96.*se_systp./systp_hat);
    CI_systp = zeros(1, 2*length(p));
    CI_systp(1:2:(2*length(p)-1)) = CI_L_systp;
    CI_systp(2:2:(2*length(p))) = CI_U_systp;
    %
    result_tp1 = [tp1_hat se_tp1 CI_tp1];
    result_tp2 = [tp2_hat se_tp2 CI_tp2];
    result_systp = [systp_hat' se_systp' CI_systp];
end



function [value_tp, value_se_tp] = system_tp_simulation(Cov_par, par, candicate_q)
    p1 = par(1); 
    a1 = par(2); b1 = par(3); 
    a2 = par(4); b2 = par(5);
    f1 = @(x) p1*wblpdf(x, a1, b1)*(1-wblcdf(x, a2, b2)) ...
             + p1*wblpdf(x, a2, b2)*(1-wblcdf(x, a1, b1)) ...
             + (1-p1)*wblpdf(x, a2, b2);
    %
    value_tp = zeros(length(candicate_q), 1);
    value_se_tp = zeros(length(candicate_q), 1);
    for i = 1:length(candicate_q)
        value_q = candicate_q(i);
        tp = tp_theta_p(par, value_q);
        %
        d_p = exp(-(tp/a2)^b2)*(1-exp(-(tp/a1)^b1));
        d_a1 = -p1*b1*a1^(-1)*(tp/a1)^b1*exp(-(tp/a1)^b1-(tp/a2)^b2);
        d_b1 = p1*(tp/a1)^b1*log(tp/a1)*exp(-(tp/a1)^b1-(tp/a2)^b2);
        d_a2 = -b2*a2^(-1)*(tp/a2)^b2*exp(-(tp/a2)^b2)*(p1*exp(-(tp/a1)^b1) + (1-p1));
        d_b2 = (tp/a2)^b2*log(tp/a2)*exp(-(tp/a2)^b2)*(p1*exp(-(tp/a1)^b1) + (1-p1));
        grad_f = [d_p d_a1 d_b1 d_a2 d_b2];
        %
        se_tp = sqrt(grad_f*Cov_par*grad_f'*(f1(tp))^(-2));
        value_se_tp(i) = se_tp;
        value_tp(i) = tp;
    end
end


function [final_t] = tp_theta_p(par, value_q)
    % CDF 
    F = @(t,x) 1-(1-x(1)).*exp(-(t./x(4)).^x(5))-x(1).*exp(-(t./x(2)).^x(3)-(t./x(4)).^x(5));
    % Grid
    if value_q<=0.1
        v1 = -8:0.00001:8;
        ds = 10.^v1;
    elseif value_q<=0.5
        v1 = -5:0.00001:8;
        ds = 10.^v1;
    else 
        v1 = -3:0.00001:8;
        ds = 10.^v1;
    end
    F_hat = F(ds,par)';
    [~,I] = min(abs(F_hat-value_q));
    grid_t = ds(I);
    final_t = grid_t;
end


