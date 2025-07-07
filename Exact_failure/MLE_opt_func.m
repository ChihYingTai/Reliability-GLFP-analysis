function [result] = MLE_opt_func(par, data)
    %
    p1 = par(1);
    par_11 = par(2); par_12 = par(3); 
    par_21 = par(4); par_22 = par(5); 
    %
    loglikelihood = 0;
    for i = 1:height(data)
        d1 = data.endtime(i);
        c1 = data.censored(i); 
        defective_idx = data.defective(i); 
        fail_idx = data.failure(i);
        if c1==0
            if ((fail_idx==1)&&(defective_idx~=2)) % 
                L1 = log( p1*wblpdf(d1, par_11, par_12)*(1-wblcdf(d1, par_21, par_22)) );
            elseif (fail_idx==2)&&(defective_idx==1)
                L1 = log( p1*wblpdf(d1, par_21, par_22)*(1-wblcdf(d1, par_11, par_12)) );
            elseif (fail_idx==2)&&(defective_idx==2)
                L1 = log( (1-p1)*wblpdf(d1, par_21, par_22) );
            elseif (fail_idx==2)&&(defective_idx==0)
                L1 = log( p1*wblpdf(d1, par_21, par_22)*(1-wblcdf(d1, par_11, par_12))...
                          + (1-p1)*wblpdf(d1, par_21, par_22) );
            else
                L1 = log(  p1*wblpdf(d1, par_11, par_12)*(1-wblcdf(d1, par_21, par_22)) ...
                         + p1*wblpdf(d1, par_21, par_22)*(1-wblcdf(d1, par_11, par_12)) ...
                         + (1-p1)*wblpdf(d1, par_21, par_22) );
            end
        else
            if defective_idx==1
                L1 = log( p1*(1-wblcdf(d1, par_11, par_12))*(1-wblcdf(d1, par_21, par_22)) );
            elseif defective_idx==2
                L1 = log( (1-p1)*(1-wblcdf(d1, par_21, par_22)) );
            else
                L1 = log( p1*(1-wblcdf(d1, par_11, par_12))*(1-wblcdf(d1, par_21, par_22)) + (1-p1)*(1-wblcdf(d1, par_21, par_22)) );
            end
        end
        loglikelihood = loglikelihood + L1;
    end
    result = -loglikelihood;



end














