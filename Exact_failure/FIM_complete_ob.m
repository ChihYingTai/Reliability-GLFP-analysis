function [Result_ob] = FIM_complete_ob(par, data)
% Observed Fisher information matrix using complete likelihood function.
% usage: [Result_ob] = FIM_complete_ob(par, data)
% 
% arguments: (input)
%  par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'endtime', 'censored', 'failure' and 'defective'.
% 
% arguments: (output)
%  Result_ob - Observed Fisher information matrix using complete likelihood function.
%
% Example usage:
%  data = readtable('...\Data\Gateoxide.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:12) = 1; data.failure(20:44) = 2;
%  data.defective = zeros(height(data), 1);
%  par = [0.4, 1.3, 0.1, 185, 8];
%  Result_complete_ob = FIM_complete_ob(par, data)
%
    %% Parameters
    p1 = par(1);
    a1 = par(2); b1 = par(3);
    a2 = par(4); b2 = par(5);
    %% Data
    data_time = data.endtime;
    %
    mi1 = (data.failure~=0)&(data.defective==0);
    mi2 = (data.failure==0)&(data.defective==0);
    mi = mi1+mi2;
    a_i = (data.failure~=0)&(data.defective~=0);
    zi = (data.failure==1);
    deltai = (data.defective==1);
    ci = (data.censored==1);
    %% Functions
    % Functions for theta
    z1_f = @(x) (x./a1).^b1;
    %
    fun_Li1 = @(x) p1.*exp(-z1_f(x)) + (1-p1);
    fun_Li2 = @(x) p1.*b1.*a1.^(-b1).*x.^(b1-1).*exp(-z1_f(x)) + b2.*a2.^(-b2).*x.^(b2-1).*fun_Li1(x);
    fun_Li3 = @(x) p1.*exp(-z1_f(x)) + (1-p1);
    %
    fun_pi2 = @(x) p1.*b1.*a1.^(-b1).*x.^(b1-1).*exp(-z1_f(x))./fun_Li2(x);
    fun_pi1 = @(x) p1.*exp(-z1_f(x))./fun_Li1(x);
    fun_pi3 = @(x) p1.*exp(-z1_f(x))./fun_Li3(x);
    %
    func_wi1 = @(x) mi2.*fun_pi2(x)                     + mi1.*zi                     + a_i.*zi;
    func_wi2 = @(x) mi2.*(1-fun_pi2(x)).*fun_pi1(x)     + mi1.*(1-zi).*fun_pi1(x)     + a_i.*(1-zi).*deltai;
    func_wi3 = @(x) mi2.*(1-fun_pi2(x)).*(1-fun_pi1(x)) + mi1.*(1-zi).*(1-fun_pi1(x)) + a_i.*(1-zi).*(1-deltai);
    func_vi1 = @(x) mi.*fun_pi3(x)     + a_i.*deltai;
    func_vi2 = @(x) mi.*(1-fun_pi3(x)) + a_i.*(1-deltai);
    %
    func_t1 = @(x) x.^(b1);
    func_t2 = @(x) x.^(b2);
    func_lt1 = @(x) log(x).*x.^(b1);
    func_lt2 = @(x) log(x).*x.^(b2);
    func_l2t1 = @(x) log(x).^(2).*x.^(b1);
    func_l2t2 = @(x) log(x).^(2).*x.^(b2);
    %% FIM
    Result_ob = zeros(5, 5); 
    Result_ob(1, 1) = p1^(-2)*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)) + ci.*func_vi1(data_time) ) ...
                      + (1-p1)^(-2)*sum( (1-ci).*func_wi3(data_time) + ci.*func_vi2(data_time) );
    %
    Result_ob(2, 2) = -b1.*a1.^(-2).*sum( (1-ci).*func_wi1(data_time) ) ...
                      + b1.*(b1+1).*a1.^(-b1-2).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_t1(data_time) ...
                                                      + ci.*func_vi1(data_time).*func_t1(data_time) );
    %
    Result_ob(3, 3) = b1.^(-2).*sum( (1-ci).*func_wi1(data_time) ) ...
                      + log(a1).^(2).*a1.^(-b1).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_t1(data_time) ...
                                                      + ci.*func_vi1(data_time).*func_t1(data_time) ) ...
                      - 2.*log(a1).*a1.^(-b1).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_lt1(data_time) ...
                                                    + ci.*func_vi1(data_time).*func_lt1(data_time) ) ...
                      + a1.^(-b1).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_l2t1(data_time) ...
                                        + ci.*func_vi1(data_time).*func_l2t1(data_time) );
    %
    Result_ob(4, 4) = - b2.*a2.^(-2).*sum( (1-ci).*(func_wi2(data_time)+func_wi3(data_time)) ) ...
                      + b2.*(b2+1).*a2.^(-b2-2).*sum( func_t2(data_time) );
    %
    Result_ob(5, 5) = b2.^(-2).*sum( (1-ci).*(func_wi2(data_time)+func_wi3(data_time)) ) ...
                      + log(a2).^(2).*a2.^(-b2).*sum( func_t2(data_time) ) ...
                      - 2.*log(a2).*a2.^(-b2).*sum( func_lt2(data_time) ) ...
                      + a2.^(-b2).*sum( func_l2t2(data_time) );
    %
    Result_ob(2, 3) = a1.^(-1).*sum( (1-ci).*func_wi1(data_time) ) ...
                      - a1.^(-b1-1).*(1-b1.*log(a1)).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_t1(data_time) ...
                                                           + ci.*func_vi1(data_time).*func_t1(data_time) ) ...
                      - b1.*a1.^(-b1-1).*sum( (1-ci).*(func_wi1(data_time)+func_wi2(data_time)).*func_lt1(data_time) ...
                                              + ci.*func_vi1(data_time).*func_lt1(data_time) );
    Result_ob(3, 2) = Result_ob(2, 3);
    %
    Result_ob(4, 5) = a2.^(-1).*sum( (1-ci).*(func_wi2(data_time)+func_wi3(data_time)) ) ...
                      - a2.^(-b2-1).*(1-b2.*log(a2)).*sum( func_t2(data_time) ) ...
                      - b2.*a2.^(-b2-1).*sum( func_lt2(data_time) );
    Result_ob(5, 4) = Result_ob(4, 5);


end









