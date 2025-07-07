function [outputArg] = FIM_missing_ob(par, data)
% Observed Fisher information matrix using missing data.
% usage: [Result_ob] = FIM_missing_ob(par, data)
% 
% arguments: (input)
%  par - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'endtime', 'censored', 'failure' and 'defective'.
% 
% arguments: (output)
%  Result_ob - Observed Fisher information matrix using missing data.
%
% Example usage:
%  data = readtable('...\Data\Gateoxide.csv');
%  data.failure = zeros(height(data), 1);
%  data.failure(1:12) = 1; data.failure(20:44) = 2;
%  data.defective = zeros(height(data), 1);
%  par = [0.4, 1.3, 0.1, 185, 8];
%  Result_missing_ob = FIM_missing_ob(par, data)
%
    %% Parameters
    p1 = par(1);
    a1 = par(2); b1 = par(3);
    a2 = par(4); b2 = par(5);
    %% Data
    data = data(((data.failure==0)|(data.defective==0)), :);
    data_time = data.endtime;
    %
    mi1 = (data.failure~=0)&(data.defective==0);
    mi2 = (data.failure==0)&(data.defective==0);
    mi = mi1+mi2;
    zi = (data.failure==1);
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
    func_t1 = @(x) x.^(b1);
    func_lt1 = @(x) log(x).*x.^(b1);
    func_l2t1 = @(x) log(x).^(2).*x.^(b1);
    %
    fun_Li2_P = @(x) b1.*a1.^(-b1).*x.^(b1-1).*exp(-z1_f(x)) + b2.*a2.^(-b2).*x.^(b2-1).*(exp(-z1_f(x)) -1);
    fun_Li2_PA1 = @(x) a1.^(-2.*b1-1).*a2.^(-b2).*b1.*x.^(b1-1).*exp(-z1_f(x)).*( b1.*a2.^(b2).*x.^(b1) + a1.^(b1).*(-b1.*a2.^(b2)+b2.*x.^(b2)) );
    fun_Li2_PB1 = @(x) a1.^(-2.*b1).*a2.^(-b2).*x.^(b1-1).*exp(-z1_f(x)).*( a1.^(b1).*a2.^(b2) + (log(a1)-log(x))...
                                                                            .*(a2.^(b2).*b1.*x.^(b1) + a1.^(b1).*(-b1.*a2.^(b2)+b2.*x.^(b2))) );
    fun_Li2_PA2 = @(x) -a2.^(-b2-1).*b2.^(2).*x.^(b2-1).*(exp(-z1_f(x))-1);
    fun_Li2_PB2 = @(x) a2.^(-b2).*x.^(b2-1).*(exp(-z1_f(x))-1).*(1-b2.*log(a2)+b2.*log(x));
    fun_Li2_A1 = @(x) p1.*a1.^(-2.*b1-1).*a2.^(-b2).*b1.*x.^(b1-1).*exp(-z1_f(x)).*( b1.*a2.^(b2).*x.^(b1) + a1.^(b1).*(-b1.*a2.^(b2)+b2.*x.^(b2)) );
    fun_Li2_A1A1 = @(x) p1.*a1.^(-3.*b1-2).*a2.^(-b2).*b1.*x.^(b1-1).*exp(-z1_f(x)).*( b1.^(2).*a2.^(b2).*x.^(2.*b1) ...
                                                                                     + a1.^(2.*b1).*(b1+1).*(b1.*a2.^(b2)-b2.*x.^(b2))...
                                                                                     - b1.*a1.^(b1).*x.^(b1).*(a2.^(b2)+3*b1.*a2.^(b2)-b2.*x.^(b2)) );
    fun_Li2_A1B1 = @(x) p1.*a1.^(-3.*b1-1).*a2.^(-b2).*x.^(b1-1).*exp(-z1_f(x)) ...
                        .*( a1.^(b1).*( 2.*a2.^(b2).*b1.*x.^(b1) + a1.^(b1).*(-2.*a2.^(b2).*b1+b2.*x.^(b2)) ) ...
                          + b1.*(log(a1)-log(x)).*( a2.^(b2).*b1.*x.^(2*b1) + a1.^(2*b1).*(a2.^(b2).*b1-b2.*x.^(b2)) ...
                                                  + a1.^(b1).*x.^(b1).*(-3.*a2.^(b2).*b1+b2.*x.^(b2)) ) );
    fun_Li2_A1A2 = @(x) -p1.*a1.^(-b1-1).*a2.^(-b2-1).*b1.*b2.^(2).*x.^(b1+b2-1).*exp(-z1_f(x));
    fun_Li2_A1B2 = @(x) p1.*a1.^(-b1-1).*a2.^(-b2).*b1.*x.^(b1+b2-1).*exp(-z1_f(x)).*(1-b2.*log(a2)+b2.*log(x));
    fun_Li2_B1 = @(x) p1.*a1.^(-2.*b1).*a2.^(-b2).*x.^(b1-1).*exp(-z1_f(x)) ...
                        .*( a1.^(b1).*a2.^(b2) + (log(a1)-log(x)).*(a2.^(b2).*b1.*x.^(b1) + a1.^(b1).*(-b1.*a2.^(b2)+b2.*x.^(b2))) );
    fun_Li2_B1B1 = @(x) p1.*a1.^(-3.*b1).*a2.^(-b2).*x.^(b1-1).*exp(-z1_f(x)).*(log(a1)-log(x)) ...
                        .*( -2.*a1.^(b1).*a2.^(b2).*(a1.^(b1)-x.^(b1)) ...
                          + (log(a1)-log(x)).*( a2.^(b2).*b1.*x.^(2*b1) + a1.^(2*b1).*(a2.^(b2).*b1-b2.*x.^(b2)) ...
                                                 + a1.^(b1).*x.^(b1).*(-3.*a2.^(b2).*b1+b2.*x.^(b2)) ) );
    fun_Li2_B1A2 = @(x) -p1.*a1.^(-b1).*a2.^(-b2-1).*b2.^(2).*x.^(b1+b2-1).*exp(-z1_f(x)).*(log(a1)-log(x));
    fun_Li2_B1B2 = @(x) p1.*a1.^(-b1).*a2.^(-b2).*x.^(b1+b2-1).*exp(-z1_f(x)).*(log(a1)-log(x)).*(1-b2.*log(a2)+b2.*log(x));
    fun_Li2_A2 = @(x) -a2.^(-b2-1).*b2.^(2).*x.^(b2-1).*( p1.*exp(-z1_f(x)) + (1-p1) );
    fun_Li2_A2A2 = @(x) a2.^(-b2-2).*b2.^(2).*(b2+1).*x.^(b2-1).*( p1.*exp(-z1_f(x)) + (1-p1) );
    fun_Li2_A2B2 = @(x) -a2.^(-b2-1).*b2.*x.^(b2-1).*( 2-b2.*log(a2)+b2.*log(x) ).*( p1.*exp(-z1_f(x)) + (1-p1) );
    fun_Li2_B2 = @(x) a2.^(-b2).*x.^(b2-1).*( p1.*exp(-z1_f(x)) + (1-p1) ).*( 1-b2.*log(a2)+b2.*log(x) );
    fun_Li2_B2B2 = @(x) -a2.^(-b2).*x.^(b2-1).*( log(a1)-log(x) ).*( p1.*exp(-z1_f(x)) + (1-p1) ).*( 2-b2.*log(a2)+b2.*log(x) );
    %
    fun_li2_PP =   @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_P(x)).*(fun_Li2_P(x)));
    fun_li2_PA1 =  @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_P(x)).*(fun_Li2_A1(x))   + (fun_Li2(x)).^(-1).*(fun_Li2_PA1(x)));
    fun_li2_PB1 =  @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_P(x)).*(fun_Li2_B1(x))   + (fun_Li2(x)).^(-1).*(fun_Li2_PB1(x)));
    fun_li2_PA2 =  @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_P(x)).*(fun_Li2_A2(x))   + (fun_Li2(x)).^(-1).*(fun_Li2_PA2(x)));
    fun_li2_PB2 =  @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_P(x)).*(fun_Li2_B2(x))   + (fun_Li2(x)).^(-1).*(fun_Li2_PB2(x)));
    fun_li2_A1A1 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A1(x)).*(fun_Li2_A1(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A1A1(x)));
    fun_li2_A1B1 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A1(x)).*(fun_Li2_B1(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A1B1(x)));
    fun_li2_A1A2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A1(x)).*(fun_Li2_A2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A1A2(x)));
    fun_li2_A1B2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A1(x)).*(fun_Li2_B2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A1B2(x)));
    fun_li2_B1B1 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_B1(x)).*(fun_Li2_B1(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_B1B1(x)));
    fun_li2_B1A2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_B1(x)).*(fun_Li2_A2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_B1A2(x)));
    fun_li2_B1B2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_B1(x)).*(fun_Li2_B2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_B1B2(x)));
    fun_li2_A2A2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A2(x)).*(fun_Li2_A2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A2A2(x)));
    fun_li2_A2B2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_A2(x)).*(fun_Li2_B2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_A2B2(x)));
    fun_li2_B2B2 = @(x) (-(fun_Li2(x)).^(-2).*(fun_Li2_B2(x)).*(fun_Li2_B2(x))  + (fun_Li2(x)).^(-1).*(fun_Li2_B2B2(x)));
    % Li3
    fun_Li3_P = @(x) exp(-z1_f(x)) - 1;
    fun_Li3_PA1 = @(x) a1.^(-b1-1).*b1.*x.^(b1).*exp(-z1_f(x));
    fun_Li3_PB1 = @(x) a1.^(-b1).*x.^(b1).*exp(-z1_f(x)).*(log(a1)-log(x));
    fun_Li3_A1 = @(x) p1.*a1.^(-b1-1).*b1.*x.^(b1).*exp(-z1_f(x));
    fun_Li3_A1A1 = @(x) -p1.*a1.^(-2.*b1-2).*b1.*x.^(b1).*exp(-z1_f(x)).*( a1.^(b1).*(b1+1)-b1.*x.^(b1) );
    fun_Li3_A1B1 = @(x) p1.*a1.^(-2.*b1-1).*x.^(b1).*exp(-z1_f(x)).*( a1.^(b1) - b1.*(a1.^(b1)-x.^(b1)).*(log(a1)-log(x)) );
    fun_Li3_B1 = @(x) p1.*a1.^(-b1).*x.^(b1).*exp(-z1_f(x)).*(log(a1)-log(x));
    fun_Li3_B1B1 = @(x) -p1.*a1.^(-2.*b1).*x.^(b1).*exp(-z1_f(x)).*(a1.^(b1)-x.^(b1)).*(log(a1)-log(x)).^(2);
    %
    fun_li3_PP =   @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_P(x)).*(fun_Li3_P(x));
    fun_li3_PA1 =  @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_P(x)).*(fun_Li3_A1(x))   + (fun_Li3(x)).^(-1).*fun_Li3_PA1(x);
    fun_li3_PB1 =  @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_P(x)).*(fun_Li3_B1(x))   + (fun_Li3(x)).^(-1).*fun_Li3_PB1(x);
    fun_li3_A1A1 = @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_A1(x)).*(fun_Li3_A1(x)) + (fun_Li3(x)).^(-1).*fun_Li3_A1A1(x);
    fun_li3_A1B1 = @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_A1(x)).*(fun_Li3_B1(x)) + (fun_Li3(x)).^(-1).*fun_Li3_A1B1(x);
    fun_li3_B1B1 = @(x) -(fun_Li3(x)).^(-2).*(fun_Li3_B1(x)).*(fun_Li3_B1(x)) + (fun_Li3(x)).^(-1).*fun_Li3_B1B1(x);
    %% FIM
    Result_ob = zeros(5, 5); 
    Result_ob(1, 1) = p1^(-2)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))...
                                 + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time) ...
                                 + mi.*ci.*fun_pi3(data_time) ) ...
                      + (1-p1)^(-2)*sum( mi2.*(1-ci).*(1-fun_pi1(data_time)).*(1-fun_pi2(data_time)) ...
                                       + mi1.*(1-ci).*(1-zi).*(1-fun_pi1(data_time)) ...
                                       + mi.*ci.*(1-fun_pi3(data_time)) ) ...
                      + sum( mi2.*(1-ci).*fun_li2_PP(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_PP(data_time) + mi.*ci.*fun_li3_PP(data_time) );
    %
    Result_ob(2, 2) = -b1*a1.^(-2)*sum( mi2.*(1-ci).*fun_pi2(data_time) ) ...
                      + b1*(b1+1)*a1^(-b1-2)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                             .*func_t1(data_time) ...
                                                + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_t1(data_time) ...
                                                + mi.*ci.*fun_pi3(data_time).*func_t1(data_time) ) ...
                      + sum( mi2.*(1-ci).*fun_li2_A1A1(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_A1A1(data_time) + mi.*ci.*fun_li3_A1A1(data_time) );
    %
    Result_ob(3, 3) = b1^(-2)*sum( mi2.*(1-ci).*fun_pi2(data_time) ) ...
                      + (log(a1))^(2)*a1^(-b1)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                               .*func_t1(data_time) ...
                                                  + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_t1(data_time) ...
                                                  + mi.*ci.*fun_pi3(data_time).*func_t1(data_time) ) ...
                      - 2*log(a1)*a1^(-b1)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                           .*func_lt1(data_time) ...
                                              + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_lt1(data_time) ...
                                              + mi.*ci.*fun_pi3(data_time).*func_lt1(data_time) ) ...
                      + a1^(-b1)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                 .*func_l2t1(data_time) ...
                                    + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_l2t1(data_time) ...
                                    + mi.*ci.*fun_pi3(data_time).*func_l2t1(data_time) ) ...
                      + sum( mi2.*(1-ci).*fun_li2_B1B1(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_B1B1(data_time) + mi.*ci.*fun_li3_B1B1(data_time) );
    %
    Result_ob(4, 4) = -b2*a2.^(-2)*sum( mi2.*(1-ci).*(1-fun_pi2(data_time)) ) + sum( mi2.*(1-ci).*fun_li2_A2A2(data_time) );
    %
    Result_ob(5, 5) = b2.^(-2)*sum( mi2.*(1-ci).*(1-fun_pi2(data_time)) ) + sum( mi2.*(1-ci).*fun_li2_B2B2(data_time) );
    %
    Result_ob(1, 2) = sum( mi2.*(1-ci).*fun_li2_PA1(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_PA1(data_time) + mi.*ci.*fun_li3_PA1(data_time) );
    Result_ob(2, 1) = Result_ob(1, 2);
    %
    Result_ob(1, 3) = sum( mi2.*(1-ci).*fun_li2_PB1(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_PB1(data_time) + mi.*ci.*fun_li3_PB1(data_time) );
    Result_ob(3, 1) = Result_ob(1, 3);
    %
    Result_ob(1, 4) = sum( mi2.*(1-ci).*fun_li2_PA2(data_time) );
    Result_ob(4, 1) = Result_ob(1, 4);
    %
    Result_ob(1, 5) = sum( mi2.*(1-ci).*fun_li2_PB2(data_time) );
    Result_ob(5, 1) = Result_ob(1, 5);
    %
    Result_ob(2, 3) = a1^(-1)*sum( mi2.*(1-ci).*fun_pi2(data_time) ) ...
                      - a1^(-b1-1)*(1-b1*log(a1))*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                                  .*func_t1(data_time) ...
                                                     + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_t1(data_time) ...
                                                     + mi.*ci.*fun_pi3(data_time).*func_t1(data_time) ) ...
                      - b1*a1^(-b1-1)*sum( mi2.*(1-ci).*(fun_pi1(data_time)+fun_pi2(data_time)-fun_pi1(data_time).*fun_pi2(data_time))... 
                                                      .*func_lt1(data_time) ...
                                         + mi1.*(1-ci).*(1-zi).*fun_pi1(data_time).*func_lt1(data_time) ...
                                         + mi.*ci.*fun_pi3(data_time).*func_lt1(data_time) ) ...
                     + sum( mi2.*(1-ci).*fun_li2_A1B1(data_time) + mi1.*(1-ci).*(1-zi).*fun_li3_A1B1(data_time) + mi.*ci.*fun_li3_A1B1(data_time) );
    Result_ob(3, 2) = Result_ob(2, 3);
    %
    Result_ob(2, 4) = sum( mi2.*(1-ci).*fun_li2_A1A2(data_time) );
    Result_ob(4, 2) = Result_ob(2, 4);
    %
    Result_ob(2, 5) = sum( mi2.*(1-ci).*fun_li2_A1B2(data_time) );
    Result_ob(5, 2) = Result_ob(2, 5);
    %
    Result_ob(3, 4) = sum( mi2.*(1-ci).*fun_li2_B1A2(data_time) );
    Result_ob(4, 3) = Result_ob(3, 4);
    %
    Result_ob(3, 5) = sum( mi2.*(1-ci).*fun_li2_B1B2(data_time) );
    Result_ob(5, 3) = Result_ob(3, 5);
    %
    Result_ob(4, 5) = a2.^(-1)*sum( mi2.*(1-ci).*(1-fun_pi2(data_time)) ) + sum( mi2.*(1-ci).*fun_li2_A2B2(data_time) );
    Result_ob(5, 4) = Result_ob(4, 5);
    %
    outputArg = Result_ob;

end




