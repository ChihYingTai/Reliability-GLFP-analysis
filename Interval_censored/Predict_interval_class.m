function [data,cmat_out] = Predict_interval_class(par_em, data, score_cut)
% Predict failure and defectiveness with confusion matrix.
% usage: [data, cmat_out] = Predict_interval_class(par_em, data, score_cut)
% 
% arguments: (input)
%  par_em - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  data - Table with columns 'starttime', 'endtime', 'censored', 
%         'failure_raw' and 'defective_raw'.
%         If 'failure_raw' and 'defective_raw' is not available, then 
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
%  par_em = [0.0115 583.804 0.274 40830 3.337];
%  data = readtable('...\Data\CB_data.csv');
%  score_cut = 0.5;
%  [data, cmat_out] = Predict_interval_class(par_em, data, score_cut)
%   
    % Parameters
    p1 = par_em(1); 
    a1 = par_em(2); b1 = par_em(3); 
    a2 = par_em(4); b2 = par_em(5);
    %
    func_g1 = @(x) p1.*wblpdf(x, a1, b1).*(1-wblcdf(x, a2, b2));
    func_g2 = @(x) p1.*wblpdf(x, a2, b2).*(1-wblcdf(x, a1, b1));
    %
    result_calss = zeros(height(data), 2);
    for i = 1:height(data)
        d1 = [data.starttime(i) data.endtime(i)];
        c1 = data.censored(i); 
        %
        value_g1 = integral(func_g1, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g2 = integral(func_g2, d1(1), d1(2), 'RelTol', 1e-8, 'AbsTol', 1e-13);
        value_g3 = (1-p1).*(wblcdf(d1(2), a2, b2)-wblcdf(d1(1), a2, b2));
        value_h1 = p1.*(1-wblcdf(d1(1), a1, b1)).*(1-wblcdf(d1(1), a2, b2));
        value_h2 = (1-p1).*(1-wblcdf(d1(1), a2, b2));
        if c1==0
            prob_infant = value_g1/(value_g1+value_g2+value_g3);
            prob_defective = (value_g1+value_g2)/(value_g1+value_g2+value_g3);
            result_calss(i, :) = [prob_infant prob_defective];
        else
            prob_defective = value_h1/(value_h1+value_h2);
            result_calss(i, :) = [0 prob_defective];
        end
    
    end
    result_calss(result_calss(1)>score_cut, 3) = 1; % class of infant
    result_calss(result_calss(2)>score_cut, 4) = 1; % class of defective
    %
    idx_noc = find(data.censored==0); 
    data.prob_infant(idx_noc) = result_calss(idx_noc, 1);
    data.class_infant(idx_noc) = result_calss(idx_noc, 3);  
    data.prob_defective = result_calss(:, 2);
    data.class_defective = result_calss(:, 4);
    %
    if any(strcmp(data.Properties.VariableNames, 'failure_raw'))&&any(strcmp(data.Properties.VariableNames, 'defective_raw'))
        cmat_out.cmat_failure = confusionmat(data.failure_raw(idx_noc)-1, 1-class_infant);
        cmat_out.cmat_defective = confusionmat(data.defective_raw-1, 1-class_defective);
    else
        cmat_out.cmat_failure = []; cmat_out.cmat_defective = [];
    end


end






