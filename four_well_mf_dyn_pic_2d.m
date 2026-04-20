clc, clear, close all;
% Parameter settings
M = 5000;               % Number of particles
mu = 5;                 % Van der Pol nonlinear parameter
sigma1 = 0.9;           % Noise intensity for x
sigma2 = 0.9;           % Noise intensity for y
T = 100;                % Total simulation time
dt = 0.01;              % Time step
N = T / dt;             % Number of time steps
frame_interval = 50;    % Record distribution every this many steps
alpha = 1;
beta = 1;

% Initial conditions
x = randn(M, 1)*0.2-3;        % Initial x values
y = randn(M, 1)*0.2-3;        % Initial y values

% Create video writer object
filename = 'particle_distribution_motion.mp4'; % Output video filename
v = VideoWriter(filename, 'MPEG-4');
v.FrameRate = 30; % Frames per second
open(v);

% Create figure window
figure;

% Time loop
for n = 1:N
    % Compute mean-field terms
    x_mean = mean(x);
    grad_W_mean_x = - (x - x_mean);
    y_mean = mean(y);           % Mean of all particles' x
    grad_W_mean_y = - (y - y_mean); % Mean-field term -2(x - mean)

    % Random increments
    dW1 = sqrt(dt) * randn(M, 1);
    dW2 = sqrt(dt) * randn(M, 1);

    % Update dynamics
    x = x - (beta * (4*x.^3 - 4*x) - alpha * grad_W_mean_x) * dt + sigma1 * dW1;
    y = y - (beta * (y.^3 - y) - alpha * grad_W_mean_y) * dt + sigma2 * dW2;
    % x = x - (beta * (x.^3 - x + 3*x.^2.*y - y + 3*x.*y.^2 + y.^3) - alpha * grad_W_mean_x) * dt + sigma1 * dW1;
    % y = y - (beta * (x.^3 - x + 3*x.^2.*y + 3*x.*y.^2 + 5*y.^3 - 5*y) - alpha * grad_W_mean_y) * dt + sigma2 * dW2;
    % x = x - (beta * (17/16*x.^3 + 45/16*x.^2.*y + 51/16*x.*y.^2+15/16*y.^3 - 5/2*x - 3/2*y) - alpha * grad_W_mean_x) * dt + sigma1 * dW1;
    % y = y - (beta * (15/16*x.^3 + 51/16*x.^2.*y + 45/16*x.*y.^2 + 17/16*y.^3 - 3/2*x - 5/2*y) - alpha * grad_W_mean_y) * dt + sigma2 * dW2;
    % x = x - (beta * (13/16 * x.^3 + (45*sqrt(3))/64 * x.^2 .* y + 27/16 * x .* y.^2 + (3*sqrt(3))/64 * y.^3 - 5/2 * x - (3*sqrt(3))/4 * y) - alpha * grad_W_mean_x) * dt + sigma1 * dW1;
    % y = y - (beta * ((15*sqrt(3))/64 * x.^3 + 27/16 * x.^2 .* y + (9*sqrt(3))/64 * x .* y.^2 + 1/16 * y.^3 - (3*sqrt(3))/4 * x - 1/2 * y) - alpha * grad_W_mean_y) * dt + sigma2 * dW2;

    % Compute and display center of mass in real time
    center_x = mean(x);
    center_y = mean(y);
    
    % Print center of mass to command line in real time
    if mod(n, frame_interval) == 0
        fprintf('Time: %.2f s, Center: (%.4f, %.4f)\n', n*dt, center_x, center_y);
    end

    % Save distribution plot every frame_interval steps
    if mod(n, frame_interval) == 0
        % Clear current figure
        clf;
        
        % Plot particle scatter
        subplot(1, 2, 1);
        scatter(x, y, 1, 'filled');
        hold on;
        % Mark center of mass on scatter plot
        plot(center_x, center_y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        hold off;
        xlabel('x'); ylabel('y');
        title(['Particle distribution - Time: ', num2str(n*dt), 's, Center: (', num2str(center_x, '%.2f'), ', ', num2str(center_y, '%.2f'), ')']);
        xlim([-5, 5]); % Adjust display range based on simulation results
        ylim([-5, 5]);
        grid on;

        % Plot 2D histogram
        subplot(1, 2, 2);
        histogram2(x, y, 30, 'FaceColor', 'flat');
        hold on;
        % Mark center of mass on histogram
        plot(center_x, center_y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        hold off;
        xlabel('x'); ylabel('y');
        title('Particle distribution - 2D histogram');
        xlim([-6, 6]); % Adjust display range based on simulation results
        ylim([-6, 6]);
        grid on;
        
        % Update figure
        drawnow;
        
        % Write current frame to video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

% Close video file
close(v);

% Completion message
disp(['Video saved as ', filename]);