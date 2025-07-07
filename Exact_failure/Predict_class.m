function [data,cmat_out] = Predict_class(par_em, data, score_cut)
% Predict failure and defectiveness with confusion matrix.
% usage: [data,cmat_out] = predict_class(par_em,data,score_cut)
% 
% arguments: (input)
%  par_em - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'endtime', 'censored', 'failure_raw' and 
%         'defective_raw'.
%         if 'failure_raw' and 'defective_raw' is not available, then 
%         cmat_out=NaN.
%  score_cut - Cut-off probability.
%         
% arguments: (output)
%  data - The data table contains additional columns: 'prob_infant',
%         'calss_infant', 'prob_defective' and 'calss_defective'.
%         'prob_infant': the probability of a unit is fail due to infant
%                        mortality.
%         'calss_infant': 1: infant mortality; 0: wear-out.
%         'prob_defective': the probability of a unit is defective.
%         'calss_defective': 1: defective; 0: non-defective.
%
% Example usage:
%  par_em = [0.4, 1.3, 0.1, 185, 8];
%  info.n_sample = 100; 
%  info.censored_t = 185.5; 
%  info.mask_candi = 0;
%  info.par = [0.1, 1.5, 0.45, 183, 7];  
%  data = Data_func(info);
%  score_cut = 0.5;
%  [data,cmat_out] = Predict_class(par_em, data, score_cut)
%   
    % Parameters
    B_p1 = par_em(1); 
    B_a1 = par_em(2); B_b1 = par_em(3); 
    B_a2 = par_em(4); B_b2 = par_em(5);
    % 
    func_g1_B = @(x) B_p1.*wblpdf(x, B_a1, B_b1).*(1-wblcdf(x, B_a2, B_b2));
    func_g2_B = @(x) B_p1.*wblpdf(x, B_a2, B_b2).*(1-wblcdf(x, B_a1, B_b1));
    func_g3_B = @(x) (1-B_p1).*wblpdf(x, B_a2, B_b2);
    func_h1_B = @(x) B_p1.*(1-wblcdf(x, B_a1, B_b1)).*(1-wblcdf(x, B_a2, B_b2));
    func_h2_B = @(x) (1-B_p1).*(1-wblcdf(x, B_a2, B_b2));
    %
    fun_pi21 = @(x) func_g1_B(x) ./ ( func_g1_B(x)+func_g2_B(x)+func_g3_B(x) );
    fun_pi22 = @(x) func_g2_B(x) ./ ( func_g1_B(x)+func_g2_B(x)+func_g3_B(x) );
    fun_pi3 = @(x) func_h1_B(x) ./ ( func_h1_B(x)+func_h2_B(x) );
    %
    [~, s1] = sort(data.endtime);
    data = data(s1, :);
    d1 = data.endtime;
    idx_noc = find(data.censored==0); idx_cen = find(data.censored==1);
    d1_failure = fun_pi21(d1(idx_noc));
    d1_defective_noc = fun_pi21(d1(idx_noc)) + fun_pi22(d1(idx_noc));
    d1_defective_cen = fun_pi3(d1(idx_cen));
    % Confusion matrix
    data.prob_defective(idx_noc) = d1_defective_noc;
    data.prob_defective(idx_cen) = d1_defective_cen;
    class_defective = ones(length(d1), 1); 
    class_defective(data.prob_defective<score_cut) = 0;
    data.class_defective = class_defective;
    %
    class_infant = ones(length(idx_noc), 1);
    class_infant(d1_failure(idx_noc)<score_cut) = 0;
    data.prob_infant(idx_noc) = d1_failure;
    data.class_infant(idx_noc) = class_infant;    
    %
    if any(strcmp(data.Properties.VariableNames, 'failure_raw'))&&any(strcmp(data.Properties.VariableNames, 'defective_raw'))
        cmat_out.cmat_failure = confusionmat(data.failure_raw(idx_noc)-1, 1-class_infant);
        cmat_out.cmat_defective = confusionmat(data.defective_raw-1, 1-class_defective);
    else
        cmat_out.cmat_failure = []; cmat_out.cmat_defective = [];
    end


end









