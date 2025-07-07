function [final_value] = EM_opt_func(initial_par, data, par_r)
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
    func_g3_B = @(x) (1-B_p1).*wblpdf(x, B_a2, B_b2);
    func_h1_B = @(x) B_p1.*(1-wblcdf(x, B_a1, B_b1)).*(1-wblcdf(x, B_a2, B_b2));
    func_h2_B = @(x) (1-B_p1).*(1-wblcdf(x, B_a2, B_b2));
    %
    fun_pi1 = @(x) func_g2_B(x) ./ ( func_g2_B(x)+func_g3_B(x) );
    fun_pi2 = @(x) func_g1_B(x) ./ ( func_g1_B(x)+func_g2_B(x)+func_g3_B(x) );
    fun_pi3 = @(x) func_h1_B(x) ./ ( func_h1_B(x)+func_h2_B(x) );
    % Functions for theta
    func_ln_g1 = @(x) log(p1) + log(wblpdf(x, a1, b1)) + log(1-wblcdf(x, a2, b2));
    func_ln_g2 = @(x) log(p1) + log(wblpdf(x, a2, b2)) + log(1-wblcdf(x, a1, b1));
    func_ln_g3 = @(x) log(1-p1) + log(wblpdf(x, a2, b2));
    func_ln_h1 = @(x) log(p1) + log(1-wblcdf(x, a1, b1)) + log(1-wblcdf(x, a2, b2));
    func_ln_h2 = @(x) log(1-p1) + log(1-wblcdf(x, a2, b2));
    %
    func_m2_lng1 = @(x) fun_pi2(x).*func_ln_g1(x);
    func_m2_lng2 = @(x) (1-fun_pi2(x)).*fun_pi1(x).*func_ln_g2(x);
    func_m2_lng3 = @(x) (1-fun_pi2(x)).*(1-fun_pi1(x)).*func_ln_g3(x);
    %
    func_m1_lng2 = @(x) fun_pi1(x).*func_ln_g2(x);
    func_m1_lng3 = @(x) (1-fun_pi1(x)).*func_ln_g3(x);
    %
    func_m_cen1 = @(x) fun_pi3(x).*func_ln_h1(x);
    func_m_cen2 = @(x) (1-fun_pi3(x)).*func_ln_h2(x);
    %% Q value
    data_time = data.endtime; 
    % No mask noc
    idx_a_noc = (data.censored==0)&(data.failure~=0)&(data.defective~=0);
    data_a_f = data_time(idx_a_noc);
    fail_a_f = 2-data.failure(idx_a_noc); 
    defect_a_f = 2-data.defective(idx_a_noc);
    value_a_noc = sum( fail_a_f.*func_ln_g1(data_a_f) + (1-fail_a_f).*defect_a_f.*func_ln_g2(data_a_f) ...
                                                      + (1-fail_a_f).*(1-defect_a_f).*func_ln_g3(data_a_f) );
    % No mask cen
    idx_a_cen = (data.censored==1)&(data.defective~=0);
    data_a_c = data_time(idx_a_cen);
    defect_a_c = 2-data.defective(idx_a_cen);
    value_a_cen = sum( defect_a_c.*func_ln_h1(data_a_c) + (1-defect_a_c).*func_ln_h2(data_a_c));
    % Mask 2 noc
    idx_m2_noc = (data.censored==0)&(data.failure==0)&(data.defective==0);
    data_m2_f = data_time(idx_m2_noc);
    value_m2_noc = sum(func_m2_lng1(data_m2_f) + func_m2_lng2(data_m2_f) + func_m2_lng3(data_m2_f));
    % Mask 1 noc
    idx_m1_noc = (data.censored==0)&(data.failure~=0)&(data.defective==0);
    data_m1_f = data_time(idx_m1_noc);
    fail_m1_f = 2-data.failure(idx_m1_noc);
    value_m1_noc = sum(fail_m1_f.*func_ln_g1(data_m1_f) + (1-fail_m1_f).*func_m1_lng2(data_m1_f) + (1-fail_m1_f).*func_m1_lng3(data_m1_f));
    % Mask cen
    idx_m_cen = (data.censored==1)&(data.defective==0);
    data_m_c = data_time(idx_m_cen);
    value_m_cen = sum(func_m_cen1(data_m_c) + func_m_cen2(data_m_c));
    %
    final_value = value_a_noc + value_a_cen + value_m2_noc + value_m1_noc + value_m_cen;
    final_value = -1*final_value;


end

