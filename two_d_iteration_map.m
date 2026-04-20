% 清除环境
clear; clc; close all;

% 设置全局字体为Times New Roman
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');

% 参数设置
kappa = 1;
% V = @(x) (x < 0) .* (x.^4/4 - x.^2/2) + (x >= 0) .* ((1.2*x).^4/4 - (1.2*x).^2/2);
V = @(x) x.^4/4 - x.^2/2 - x.^3/3;

% lambda 和 m 的范围
lambda_min = -40;
lambda_max = -0.05;
m_min = -2;
m_max = 2;

% 创建网格
lambda_vals = linspace(lambda_min, lambda_max, 200);
m_vals = linspace(m_min, m_max, 100);
[LAMBDA, M] = meshgrid(lambda_vals, m_vals);

% 预分配数组
F1_vals = zeros(size(LAMBDA));
F2_vals = zeros(size(LAMBDA));

% 积分范围
integral_range = 5;

% 循环计算每个 (lambda, m) 对应的 F1 和 F2
total_points = numel(LAMBDA);
current_point = 0;

fprintf('开始计算二维映射 F(lambda, m)...\n');

for i = 1:size(LAMBDA, 1)
    for j = 1:size(LAMBDA, 2)
        lambda = LAMBDA(i, j);
        m = M(i, j);
        current_point = current_point + 1;
        
        try
            % 计算分母的积分 (配分函数)
            denom_integral = integral(@(x) exp(lambda * (V(x)/kappa + x.^2/2 - m * x)), ...
                                     -integral_range, integral_range, 'RelTol', 1e-6, 'AbsTol', 1e-8);
            
            % 计算一阶矩的积分
            num_integral_F1 = integral(@(x) x .* exp(lambda * (V(x)/kappa + x.^2/2 - m * x)), ...
                                      -integral_range, integral_range, 'RelTol', 1e-6, 'AbsTol', 1e-8);
            
            % 计算二阶矩的积分
            num_integral_F2 = integral(@(x) x.^2 .* exp(lambda * (V(x)/kappa + x.^2/2 - m * x)), ...
                                      -integral_range, integral_range, 'RelTol', 1e-6, 'AbsTol', 1e-8);
            
            % 计算 F1 (均值)
            if abs(denom_integral) > 1e-12
                F1_vals(i, j) = num_integral_F1 / denom_integral;
                
                % 计算 F2 (负的方差的倒数)
                variance = (num_integral_F2 / denom_integral) - (F1_vals(i, j))^2;
                if abs(variance) > 1e-12
                    F2_vals(i, j) = -1 / variance;
                else
                    F2_vals(i, j) = NaN;
                end
            else
                F1_vals(i, j) = NaN;
                F2_vals(i, j) = NaN;
            end
            
        catch ME
            % 如果积分失败，使用 NaN
            F1_vals(i, j) = NaN;
            F2_vals(i, j) = NaN;
        end
        
        % 显示进度
        if mod(current_point, 100) == 0
            fprintf('已完成 %.1f%%\n', current_point/total_points*100);
        end
    end
end

fprintf('计算完成！\n');

% 创建图形
figure('Position', [100, 100, 600, 500]);

% 计算 F1(lambda, m) = m 和 F2(lambda, m) = lambda 的交线
F1_minus_m = F1_vals - M;
F2_minus_lambda = F2_vals - LAMBDA;

% 提取 F1 = m 的等高线
C1 = contourc(lambda_vals, m_vals, F1_minus_m, [0 0]);
% 提取 F2 = lambda 的等高线
C2 = contourc(lambda_vals, m_vals, F2_minus_lambda, [0 0]);

% 绘制在定义域上的两条交线投影
hold on;

% 初始化曲线绘制句柄
h_blue = [];
h_red = [];

% 绘制 F1 = m 的投影 (蓝色曲线)
if ~isempty(C1)
    idx = 1;
    while idx < size(C1, 2)
        level = C1(1, idx);
        n_points = C1(2, idx);
        if level == 0 && n_points > 1
            contour_x1 = C1(1, idx+1:idx+n_points);
            contour_y1 = C1(2, idx+1:idx+n_points);
            h_blue = plot(contour_x1, contour_y1, 'k-', 'LineWidth', 2, 'DisplayName', 'blue');
            
            % 存储第一条曲线的所有点
            curve1_points = [contour_x1' contour_y1'];
        end
        idx = idx + n_points + 1;
    end
end

% 绘制 F2 = lambda 的投影 (红色曲线)
if ~isempty(C2)
    idx = 1;
    while idx < size(C2, 2)
        level = C2(1, idx);
        n_points = C2(2, idx);
        if level == 0 && n_points > 1
            contour_x2 = C2(1, idx+1:idx+n_points);
            contour_y2 = C2(2, idx+1:idx+n_points);
            h_red = plot(contour_x2, contour_y2, 'r-.', 'LineWidth', 2, 'DisplayName', 'red');
            
            % 存储第二条曲线的所有点
            curve2_points = [contour_x2' contour_y2'];
        end
        idx = idx + n_points + 1;
    end
end

% 设置坐标轴标签，纵轴使用LaTeX格式
xlabel('\lambda', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('$m$', 'FontSize', 14, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

% 设置标题
title('Bifurcation Structure for the Asymmetric Double-well Model', ...
      'FontSize', 16, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

% 设置图例 - 显示红蓝两个图例
% 创建图例句柄数组
legend_handles = [];
legend_labels = {};

if ~isempty(h_blue)
    legend_handles = [legend_handles, h_blue];
    legend_labels = [legend_labels, {'Solutions for the moment bifurcation equation'}];
end

if ~isempty(h_red)
    legend_handles = [legend_handles, h_red];
    legend_labels = [legend_labels, {"Solutions for Dawson's criterion"}];
end

% 显示图例
if ~isempty(legend_handles)
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 12, 'FontName', 'Times New Roman');
end

% 设置网格和坐标范围
grid on;
xlim([lambda_min, lambda_max]);
ylim([m_min, m_max]);

% 设置坐标轴刻度字体
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

fprintf('图像绘制完成！\n');