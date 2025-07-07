function Hazard_plot(par_em, plot_info)
% Generate hazard plot.
% usage: hazard_plot_gate(par_em, plot_info)
% 
% arguments: (input)
%  par_em - vector of (\pi,\alpha_1,\beta_1,\alpha_2,\beta_2).
%  plot_info - Structure with transform: 'nature_log' for gate oxide data 
%                                         and 'log_10' for CB data;
%                             points: plot points;
%                             x_label: labels of x-axis of log transformation;
%                             y_label: labels of y-axis of log transformation;
%                             x_lim: the x-axis limits;
%                             y_lim: the y-axis limits.
%         
% arguments: (output)
%  hazaed plot.
%
% Example usage:
%  % Exact failure
%  par_em = [0.483 3.457 0.118 185.939 8.407];
%  plot_info.transform = 'nature_log'
%  plot_info.points = min(data.endtime):0.01:250;
%  plot_info.x_label = -6:1.5:6;
%  plot_info.y_label = -8:1:1;
%  plot_info.x_lim = [-6 6];
%  plot_info.y_lim = [-8 1];
%  Hazard_plot(par_em, plot_info)
%
%  % Interval data
%  plot_info.transform = 'log_10'
%  plot_info.points = 0:10^5;
%  plot_info.x_label = 0:10;
%  plot_info.y_label = -8:-3;
%  plot_info.x_lim = [0 5.2];
%  plot_info.y_lim = [-8 -3];
%
%   
    %% Settings
    B_p1 = par_em(1); 
    B_a1 = par_em(2); B_b1 = par_em(3); 
    B_a2 = par_em(4); B_b2 = par_em(5);
    %
    func_g1_B = @(x) B_p1.*wblpdf(x, B_a1, B_b1).*(1-wblcdf(x, B_a2, B_b2));
    func_g2_B = @(x) B_p1.*wblpdf(x, B_a2, B_b2).*(1-wblcdf(x, B_a1, B_b1));
    func_g3_B = @(x) (1-B_p1).*wblpdf(x, B_a2, B_b2);
    func_h1_B = @(x) B_p1.*(1-wblcdf(x, B_a1, B_b1)).*(1-wblcdf(x, B_a2, B_b2));
    func_h2_B = @(x) (1-B_p1).*(1-wblcdf(x, B_a2, B_b2));
    fun_hi = @(x) (func_g1_B(x)+func_g2_B(x)+func_g3_B(x)) ./ ( func_h1_B(x)+func_h2_B(x) );
    %% Hazard plot  
    if plot_info.transform=='nature_log'
        d1 = plot_info.points;
        d1 = d1(2:end);
        x1 = log(d1);
        plot_d1 = log(fun_hi(d1));
        %
        t = tiledlayout(1,1);
        ax1 = axes(t);
        hold on
        plot(x1, plot_d1, 'k-','LineWidth', 1.2)
        x_label = plot_info.x_label; y_label = plot_info.y_label;
        set(ax1, 'XTick', x_label, 'XTickLabel', string(round(exp(x_label), 2)), ...
                 'YTick', y_label, 'YTickLabel', string(round(exp(y_label), 3)), ...
                 'XLim', plot_info.x_lim, 'YLim', plot_info.y_lim, ...
                 'Box', 'on');
        ax1.XLabel.String = 'Time (Seconds)';
        ax1.YLabel.String = 'Hazard';
    elseif plot_info.transform=='log_10'
        d1 = plot_info.points;
        d1 = d1(2:end);
        x1 = log10(d1);
        plot_d1 = log10(fun_hi(d1))';
        %
        t = tiledlayout(1,1);
        ax1 = axes(t);
        hold on
        plot(x1, plot_d1, 'k-','LineWidth', 1.2)
        x_label = plot_info.x_label; 
        x_tick_labels = arrayfun(@(x) ['10^{', num2str(x), '}'], x_label, 'UniformOutput', false);
        x_tick_labels = [x_tick_labels, {'interpreter', 'latex'}];
        y_label = plot_info.y_label;
        y_tick_labels = arrayfun(@(x) ['10^{', num2str(x), '}'], y_label, 'UniformOutput', false);
        y_tick_labels = [y_tick_labels, {'interpreter', 'latex'}];
        set(ax1, 'XTick', x_label, 'XTickLabel', x_tick_labels, ...
                 'YTick', y_label, 'YTickLabel', y_tick_labels, ...
                 'XLim', plot_info.x_lim, 'YLim', plot_info.y_lim, ...
                 'Box', 'on');
        ax1.XLabel.String = 'Time (Hours)';
        ax1.YLabel.String = 'Hazard';
    end

end







