function Probability_plot(par_em, cov, data, plot_info)
% Generate probability plot.
% usage: probability_plot_gate(par_em, cov, data)
% 
% arguments: (input)
%  par_em - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  cov - covariance matrix of par_em.
%  data - Table with columns 'endtime', 'censored'.  
%  plot_info - Structure with points: plot points;
%                             x_label: labels of x-axis of log transformation;
%                             y_tex: probability labels in y-axis;
%                             x_lim: the x-axis limits;
%                             y_lim: the y-axis limits.
% 
% arguments: (output)
%  probability plot of gate oxide data.
%
% Example usage:
%  % Exact failure
%  par_em = [0.4, 1.3, 0.1, 185, 8];
%  FIM_c = FIM_complete_ob(par_em, data);
%  FIM_m = FIM_missing_ob(par_em, data);
%  cov = inv(FIM_c - FIM_m);
%  v1 = -11:0.001:3;
%  plot_info.points = 10.^v1;
%  plot_info.x_label = -10:2:2;
%  plot_info.y_tex = [0.001,0.003,0.005,0.01,0.02,0.03,0.04,0.06,0.1,0.15,0.2,0.3,0.4,0.5,0.7,0.8,0.9,0.96,0.99];
%  plot_info.x_lim = [-10 3];
%  plot_info.y_lim = [-2.5 0.5];
%  Probability_plot(par_em, cov, data, plot_info)
%
%  % Interval data
%  v1 = 0:0.001:6;
%  plot_info.points = 10.^v1;
%  plot_info.x_label = 0:6;
%  plot_info.y_tex = [0.001,0.003,0.005,0.01,0.02,0.03,0.06,0.1,0.2,0.3,0.5,0.7,0.9,0.96];
%  plot_info.x_lim = [-0.2 4.65];
%  plot_info.y_lim = [-3.5 0.25];
%   
    %% Settings
    p1 = par_em(1); 
    a1 = par_em(2); b1 = par_em(3); 
    a2 = par_em(4); b2 = par_em(5);
    F = @(x) 1-(1-p1).*exp(-(x./a2).^b2)-p1.*exp(-(x./a1).^b1-(x./a2).^b2);
    F_P1 = @(x) -exp(-(x./a2).^b2).*(exp(-(x./a1).^b1)-1);
    F_A1 = @(x) -p1.*b1.*a1.^(-b1-1).*x.^(b1).*exp(-(x./a1).^b1-(x./a2).^b2);
    F_B1 = @(x) -p1.*a1.^(-b1).*x.^(b1).*(log(a1)-log(x)).*exp(-(x./a1).^b1-(x./a2).^b2);
    F_A2 = @(x) -a2.^(-b2-1).*b2.*x.^(b2).*exp(-(x./a2).^b2).*( p1.*exp(-(x./a1).^b1)+(1-p1) );
    F_B2 = @(x) -a2.^(-b2).*x.^(b2).*(log(a2)-log(x)).*( p1.*exp(-(x./a1).^b1)+(1-p1) ).*exp(-(x./a2).^b2);
    ds = plot_info.points;
    g1 = [F_P1(ds)', F_A1(ds)', F_B1(ds)', F_A2(ds)', F_B2(ds)'];
    %% F_hat and CI - Ic-Im
    F_hat = F(ds)';
    se_F = sqrt(sum(g1*cov.*g1, 2));
    CI_L = F_hat./( F_hat+(1-F_hat).*exp(1.96.*se_F./(F_hat.*(1-F_hat))) );
    CI_U = F_hat./( F_hat+(1-F_hat).*exp(-1.96.*se_F./(F_hat.*(1-F_hat))) );
    s1 = find(diff(CI_L)<0, 1, 'first');
    CI_L(s1:end) = CI_L(s1-1);
    %% F_hat nonparametric
    failure_idx = find(data.censored==0);
    if sum(strcmp(data.Properties.VariableNames, 'count'))==0 
        F_hat_np = (1:length(failure_idx))./height(data);
    else
        F_hat_np = cumsum(data1.count)./sum(data1.count);
        F_hat_np = F_hat_np(1:end-1);
    end
    %% Probability plot
    x = log10(data.endtime(failure_idx));
    y = log10(-log10(1-F_hat_np));
    x1 = log10(ds); 
    z1 = log10(-log10(1-F_hat));
    s2 = 1:length(ds); % [1:100:13100 13060:10:13700 13700:25:length(ds)];  
    y3 = log10(-log10(1-CI_L)); y4 = log10(-log10(1-CI_U));
    % Plot first axes
    t = tiledlayout(1,1);
    ax1 = axes(t);
    plot(x, y, 'k.', 'MarkerSize', 5)
    x_label = plot_info.x_label;
    y_tex = plot_info.y_tex;
    y_label = log10(-log10(1-y_tex));
    % string(x_label)
    x_tick_labels = arrayfun(@(x) ['10^{', num2str(x), '}'], x_label, 'UniformOutput', false);
    x_tick_labels = [x_tick_labels, {'interpreter', 'latex'}];
    set(ax1, 'XTick', x_label, 'XTickLabel', x_tick_labels, ... 
             'YTick', y_label, 'YTickLabel', string(y_tex), ...    
             'XLim', plot_info.x_lim, 'YLim', plot_info.y_lim, ... 
             'Box', 'on');%
    ax1.XLabel.String = 'Seconds';
    ax1.YLabel.String = 'Proportion Failing';
    % Plot second axes
    ax2 = axes(t);
    hold on
    plot(x, y, 'k.', 'MarkerSize', 20)
    plot(x1, z1, 'k.', 'MarkerSize', 1)
    plot(x1(s2), y3(s2), 'k.', 'MarkerSize', 2)
    plot(x1(s2), y4(s2), 'k.', 'MarkerSize', 2)
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    set(ax2, 'XLim', plot_info.x_lim, 'YLim', plot_info.y_lim, ... 
             'Box', 'on'); 
    ax2.XLabel.String = 'Log Seconds';
    ax2.YLabel.String = 'Standard Quantile';

end











