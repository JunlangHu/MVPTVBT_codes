clear; clc; close all;

%% Parameter settings
t_range = linspace(0, 1, 101);  % homotopy parameter t from 0 to 1
m0 = 0;                         % initial m (symmetric case)
lambda0 = -2.2;                 % initial lambda (theoretical value)

% Store solution trajectories
m_traj = zeros(size(t_range));
lambda_traj = zeros(size(t_range));
m_traj(1) = m0;
lambda_traj(1) = lambda0;

%% Numerical integration settings
x_min = -3; x_max = 3; nx = 2001;  % increased points to better detect extrema
x = linspace(x_min, x_max, nx);
dx = x(2) - x(1);

%% Fixed planar region settings
m_limits = [-1, 1];            % fixed range for m
lambda_limits = [-20, -1];     % fixed range for lambda

%% Homotopy continuation - predictor-corrector method
fprintf('Starting homotopy continuation...\n');

for i = 1:length(t_range)-1
    t = t_range(i);
    m_current = m_traj(i);
    lambda_current = lambda_traj(i);
    
    % Predictor step (Euler method)
    [H, J, ~, ~] = compute_H_J(m_current, lambda_current, t, x, dx);
    H_t = compute_H_t(m_current, lambda_current, t, x, dx);
    
    if rcond(J) < 1e-12
        fprintf('Jacobian singular at t=%.3f, stopping computation\n', t);
        break;
    end
    
    step = t_range(i+1) - t_range(i);
    predictor = -J \ H_t;
    m_pred = m_current + step * predictor(1);
    lambda_pred = lambda_current + step * predictor(2);
    
    % Corrector step (Newton's method)
    max_newton_iter = 10;
    tol = 1e-8;
    
    for iter = 1:max_newton_iter
        [H_pred, J_pred, E_x, Var_x] = compute_H_J(m_pred, lambda_pred, t_range(i+1), x, dx);
        
        if norm(H_pred) < tol
            break;
        end
        
        if rcond(J_pred) < 1e-12
            fprintf('Jacobian singular in correction step\n');
            break;
        end
        
        correction = -J_pred \ H_pred;
        m_pred = m_pred + correction(1);
        lambda_pred = lambda_pred + correction(2);
    end
    
    m_traj(i+1) = m_pred;
    lambda_traj(i+1) = lambda_pred;
    
    fprintf('t=%.3f: m=%.6f, lambda=%.6f, |H|=%.2e\n', ...
            t_range(i+1), m_pred, lambda_pred, norm(H_pred));
end

%% Improved extrema detection function
function [max_peaks, max_locs, min_peaks, min_locs] = find_all_extrema(p, x, min_prominence)
    % Find all local maxima
    [max_peaks, max_locs] = findpeaks(p, x, 'MinPeakProminence', min_prominence);
    
    % If not enough maxima are found, lower the requirement and try again
    % if length(max_peaks) < 2
    %     [max_peaks_temp, max_locs_temp] = findpeaks(p, x, 'SortStr', 'descend');
    %     if ~isempty(max_peaks_temp)
    %         % Merge results and remove duplicates
    %         all_max_peaks = [max_peaks; max_peaks_temp];
    %         all_max_locs = [max_locs; max_locs_temp];
    %         [~, unique_idx] = uniquetol(all_max_locs, 1e-3);
    %         max_peaks = all_max_peaks(unique_idx);
    %         max_locs = all_max_locs(unique_idx);
    %     end
    % end
    
    % Sort by location
    % if ~isempty(max_peaks)
    %     [max_locs, sort_idx] = sort(max_locs);
    %     max_peaks = max_peaks(sort_idx);
    % end
    
    % Find all local minima
    [min_peaks, min_locs] = findpeaks(-p, x, 'MinPeakProminence', min_prominence/2);
    min_peaks = -min_peaks;
    
    if ~isempty(min_peaks)
        [min_locs, sort_idx] = sort(min_locs);
        min_peaks = min_peaks(sort_idx);
    end
end

%% Create animation
animation_fig = figure('Position', [100, 100, 1200, 500]);
frame_counter = 1;

for i = 1:2:length(t_range)  % display every other frame for smoother animation
    if i > length(m_traj) || isnan(m_traj(i))
        break;
    end
    
    figure(animation_fig);
    clf;
    
    % Subplot 1: Trajectory (fixed region)
    subplot(1,2,1);
    hold on;
    grid on;
    
    % Plot full trajectory
    plot(m_traj(1:i), lambda_traj(1:i), 'b-', 'LineWidth', 2);
    % Plot current point
    plot(m_traj(i), lambda_traj(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Set fixed region
    xlim(m_limits);
    ylim(lambda_limits);
    
    xlabel('m');
    ylabel('\lambda');
    title(sprintf('Fixed point trajectory (t=%.2f)', t_range(i)));
    
    % Mark start and end points
    if i == 1
        plot(m_traj(1), lambda_traj(1), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        legend('Trajectory', 'Current point', 'Start point', 'Location', 'best');
    elseif i >= length(m_traj)-2
        plot(m_traj(1), lambda_traj(1), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        plot(m_traj(end), lambda_traj(end), 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
        legend('Trajectory', 'Current point', 'Start point', 'End point', 'Location', 'best');
    else
        legend('Trajectory', 'Current point', 'Location', 'best');
    end
    
    % Subplot 2: Current probability density and extrema
    subplot(1,2,2);
    hold on;
    grid on;
    
    [~, ~, E_x, Var_x, p] = compute_H_J(m_traj(i), lambda_traj(i), t_range(i), x, dx);
    
    % Plot probability density function
    h_p = plot(x, p, 'r-', 'LineWidth', 2);
    
    % Find all extrema (lower prominence threshold to detect more points)
    min_prominence = 0.01;  % lowered prominence threshold
    [max_peaks, max_locs, min_peaks, min_locs] = find_all_extrema(p, x, min_prominence);
    
    % Mark all local maxima
    if ~isempty(max_peaks)
        h_max = plot(max_locs, max_peaks, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'green');
        
        % Display coordinates of all maxima
        for j = 1:length(max_peaks)
            % Adjust text position to avoid overlap
            if max_locs(j) < 0
                text_pos = [max_locs(j)-0.3, max_peaks(j)+0.05];
                halign = 'right';
            else
                text_pos = [max_locs(j)+0.3, max_peaks(j)+0.05];
                halign = 'left';
            end
            
            text(text_pos(1), text_pos(2), sprintf('(%.2f,%.3f)', max_locs(j), max_peaks(j)), ...
                 'FontSize', 8, 'HorizontalAlignment', halign, 'Color', 'green', ...
                 'BackgroundColor', 'white', 'EdgeColor', 'green');
        end
    else
        h_max = [];
    end
    
    % Mark all local minima
    if ~isempty(min_peaks)
        h_min = plot(min_locs, min_peaks, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'green');
        
        % Display only significant minima (depth > half of max peak)
        significant_min = min_peaks < 0.5 * max(p);
        for j = 1:length(min_peaks)
            if significant_min(j)
                text(min_locs(j), min_peaks(j)-0.05, sprintf('(%.2f,%.3f)', min_locs(j), min_peaks(j)), ...
                     'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', 'green', ...
                     'BackgroundColor', 'white', 'EdgeColor', 'green');
            end
        end
    else
        h_min = [];
    end
    
    xlabel('x');
    ylabel('p(x)');
    title('Probability density function and extrema');
    xlim([-2.5, 2.5]);
    ylim([0, 2]);
    
    % Add statistical information
    text(-2.4, 1.1, sprintf('E[x]=%.3f\nVar[x]=%.3f', E_x, Var_x), ...
         'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Dynamically create legend
    legend_items = {'p(x)'};
    if ~isempty(h_max)
        legend_items{end+1} = 'Local maxima';
    end
    if ~isempty(h_min)
        legend_items{end+1} = 'Local minima';
    end
    legend(legend_items, 'Location', 'northeast');
    
    drawnow;
    
    % Save frame for animation
    frame(frame_counter) = getframe(animation_fig);
    frame_counter = frame_counter + 1;
end

%% Save animation
fprintf('Saving animation...\n');
if exist('frame', 'var') && ~isempty(frame)
    video_writer = VideoWriter('homotopy_continuation_detailed', 'MPEG-4');
    video_writer.FrameRate = 5;
    open(video_writer);
    for i = 1:length(frame)
        if i <= length(frame) && ~isempty(frame(i).cdata)
            writeVideo(video_writer, frame(i));
        end
    end
    close(video_writer);
    fprintf('Animation saved as homotopy_continuation_detailed.mp4\n');
else
    fprintf('No valid frames to save\n');
end

%% Plot final results
figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
plot(t_range(1:length(m_traj)), m_traj, 'b-', 'LineWidth', 2);
grid on;
xlabel('Homotopy parameter t');
ylabel('m');
title('Evolution of m with t');
xlim([0, 1]);

subplot(1,2,2);
plot(t_range(1:length(lambda_traj)), lambda_traj, 'r-', 'LineWidth', 2);
grid on;
xlabel('Homotopy parameter t');
ylabel('\lambda');
title('Evolution of \lambda with t');
xlim([0, 1]);

%% Display detailed extrema evolution information
fprintf('\n=== Extremum Evolution Analysis ===\n');
sample_points = [1, round(length(t_range)/3), round(2*length(t_range)/3), length(t_range)];

for idx = 1:length(sample_points)
    i = sample_points(idx);
    if i > length(m_traj), continue; end
    
    fprintf('\nTime t=%.2f:\n', t_range(i));
    [~, ~, ~, ~, p] = compute_H_J(m_traj(i), lambda_traj(i), t_range(i), x, dx);
    
    [max_peaks, max_locs, min_peaks, min_locs] = find_all_extrema(p, x, 0.01);
    
    fprintf('Local maxima (%d):\n', length(max_peaks));
    for j = 1:length(max_peaks)
        fprintf('  Peak %d: x=%.4f, p=%.4f\n', j, max_locs(j), max_peaks(j));
    end
    
    fprintf('Local minima (%d):\n', length(min_peaks));
    for j = 1:length(min_peaks)
        fprintf('  Valley %d: x=%.4f, p=%.4f\n', j, min_locs(j), min_peaks(j));
    end
end

%% Helper function: compute H and Jacobian
function [H, J, E_x, Var_x, p] = compute_H_J(m, lambda, t, x, dx)
    % Compute potential energy function
    V = x.^4/4 - x.^2/2 - (t/3)*x.^3;
    V = V/1;
    
    % Compute exponent part
    exponent = lambda * (V + x.^2/2 - m*x);
    
    % Avoid numerical overflow
    max_exp = max(exponent);
    exp_weight = exp(exponent - max_exp);
    
    % Compute partition function
    Z = sum(exp_weight) * dx;
    
    % Probability density
    p = exp_weight / Z;
    
    % Expected values
    E_x = sum(x .* p) * dx;
    E_x2 = sum(x.^2 .* p) * dx;
    E_x3 = sum(x.^3 .* p) * dx;
    E_x4 = sum(x.^4 .* p) * dx;
    
    % Variance and third central moment
    Var_x = E_x2 - E_x^2;
    mu3 = E_x3 - 3*E_x*E_x2 + 2*E_x^3;
    
    % A(x) = V(x;t) + x^2/2 - m*x
    A_x = V + x.^2/2 - m*x;
    E_A = sum(A_x .* p) * dx;
    
    % Covariance calculations
    Cov_x_A = sum(x .* A_x .* p) * dx - E_x * E_A;
    Cov_x2_A = sum(x.^2 .* A_x .* p) * dx - E_x2 * E_A;
    
    % H function
    H1 = m - E_x;
    H2 = lambda + 1/Var_x;
    H = [H1; H2];
    
    % Jacobian matrix
    J11 = 1 + lambda * Var_x;
    J12 = -Cov_x_A;
    J21 = (lambda * mu3) / (Var_x^2);
    J22 = 1 - (Cov_x2_A - 2*E_x*Cov_x_A) / (Var_x^2);
    
    J = [J11, J12; J21, J22];
end

%% Helper function: compute derivative of H with respect to t
function H_t = compute_H_t(m, lambda, t, x, dx)
    % Compute current probability distribution
    [~, ~, E_x, Var_x, p] = compute_H_J(m, lambda, t, x, dx);
    
    % Compute moments
    E_x = sum(x .* p) * dx;
    E_x2 = sum(x.^2 .* p) * dx;
    E_x3 = sum(x.^3 .* p) * dx;
    E_x4 = sum(x.^4 .* p) * dx;
    
    % Covariance calculations
    Cov_x_x3 = E_x4 - E_x * E_x3;
    Cov_x2_x3 = sum(x.^5 .* p) * dx - E_x2 * E_x3;
    
    % H_t calculation
    H_t1 = (lambda/3) * Cov_x_x3;
    H_t2 = (lambda/(3*Var_x^2)) * (Cov_x2_x3 - 2*E_x*Cov_x_x3);
    
    H_t = [H_t1; H_t2];
end