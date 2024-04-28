% MATLAB script written by Salvatore Granata
% This MATLAB script contains the exercise 3 and, in another section, part
% of the exercise 4 (which require the ex3 to work itself).

close all
clear all
clc

qi = [0; 0; 0];     % initial configuration

% Random final configuration s.t. ||qf - qi|| = 1
qf = rand(3, 1);
qf = (qf - qi) ./ norm(qf - qi);

% Time parameters
Ts = 0.001;             %sample time
trajDuration = 1;       %trajectory duration
t = 0:Ts:trajDuration;  %time vector

% arrays used to store the trajectory
s = zeros(size(t));
ds = zeros(size(t));
dds = zeros(size(t));

% definition of cubic polynomial coefficients
a0 = 0;
a1 = 0;
a2 = -3 * a0 / trajDuration^2 - 2 * a1 / trajDuration + 3 / trajDuration^2 - 0 / trajDuration;
a3 = 2 * a0 / trajDuration^3 + a1 / trajDuration^2 - 2 / trajDuration^3 + 0 / trajDuration^2;

% trajectory computation with cubic polynomial function (implemented below)
for i = 1:length(t)
    [s(i), ds(i), dds(i)] = cubic_polinomial(a0, a1, a2, a3, t(i));
end

k = 1;           
theta_i = qi(3);
theta_f = qf(3);
x_i = qi(1);
y_i = qi(2);
x_f = qf(1);
y_f = qf(2);

alpha_x = k * cos(theta_f) - 3 * x_f;
alpha_y = k * sin(theta_f) - 3 * y_f;
beta_x = k * cos(theta_i) + 3 * x_i;
beta_y = k * sin(theta_i) + 3 * y_i;

% positions along the path
x_s = s.^3 * x_f - (s - 1).^3 * x_i + alpha_x * s.^2 .* (s - 1) + beta_x * s .* (s - 1).^2;
y_s = s.^3 * y_f - (s - 1).^3 * y_i + alpha_y * s.^2 .* (s - 1) + beta_y * s .* (s - 1).^2;

% derivatives of positions
x_p=3*x_f*s.^2-3*x_i*(s-1).^2+alpha_x*s.^2+2*alpha_x*(s-1).*s+beta_x*(s-1).^2+2*beta_x*s.*(s-1);
y_p = 3*y_f*s.^2-3*y_i*(s-1).^2+alpha_y*s.^2+2*alpha_y*(s-1).*s+beta_y*(s-1).^2+2*beta_y*s.*(s-1);

% double derivatives of positions
x_pp=6*x_f*s-3*x_i*(s-1)+2*alpha_x*(s-1)+4*alpha_x*s+4*beta_x*(s-1)+2*beta_x*s;
y_pp=6*y_f*s-3*y_i*(s-1)+2*alpha_y*(s-1)+4*alpha_y*s+4*beta_y*(s-1)+2*beta_y*s;

% theta angle computation
theta = atan2(y_p,x_p);

v_s = sqrt(x_p.^2 + y_p.^2);                            %linear velocity v(s)
w_s = (y_pp .* x_p - x_pp .* y_p) ./ (x_p.^2 + y_p.^2); %angular velocity v(s)

v_t = v_s .* ds;    %linear velocity v(t)
w_t = w_s .* ds;    %angular velocity w(t)


% Plot of the trajectory
figure;
plot(x_s, y_s, 'LineWidth', 2);
hold on;
plot(qi(1), qi(2), 'ro', 'MarkerSize', 10); % initial configuration
plot(qf(1), qf(2), 'go', 'MarkerSize', 10); % final configuration
xlabel('x (m)');
ylabel('y (m)');
title('Robot Path');
legend('Path', 'qi', 'qf');
axis equal;
hold off;

% Plot of the cubic trajectory
figure;
subplot(3,1,1);
plot(t, s);                 %trajectory expected to be cubic
xlabel('Time (s)');
ylabel('Position (m)');
title('Position along the path');

subplot(3,1,2)
plot(t,ds);                 %velocity expected to be parabolic
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity along the path');

subplot(3,1,3);
plot(t, dds);               %acceleration expected to be linear
xlabel('Time (s)');
ylabel('Linear Acceleration (m/s^2)');
title('Acceleration along the path');

figure      
plot(t,theta);       %theta
xlabel('Times (s)');
ylabel('Theta');
title('Robot angular position');

figure;
subplot(2,1,1);
plot(t, v_s);               % plot of linear velocity (useful to check if constraint is satisfied)
xlabel('Time (s)');
ylabel('Linear Velocity (m/s)');
title('Linear Velocity v(s)');

subplot(2,1,2);
plot(t, w_s);               % plot of angular velocity (useful to check if constraint is satisfied)
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity w(s)');

figure;
subplot(2,1,1);
plot(t, v_t);               % plot of linear velocity (useful to check if constraint is satisfied)
xlabel('Time (s)');
ylabel('Linear Velocity (m/s)');
title('Linear Velocity v(t)');

subplot(2,1,2);
plot(t, w_t);               % plot of angular velocity (useful to check if constraint is satisfied)
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity w(t)');

% constraints
v_max = 2;      % maximum linear velocity
w_max = 1;      % maximum angular velocity

v_scaling = v_max / max(abs(v_t));
w_scaling = w_max / max(abs(w_t));
scaling = min(v_scaling, w_scaling);

% Check if velocity constraints are violated
if scaling < 1

    % If constraints are violated, repeat the trajectory computation with scaled time
    disp('Velocity constraints are not satisfied. Time scaling is needed.');
    t_scaled = t / scaling; %minimum time s.t. the constraints are satisfied 
                                   
    t=t_scaled;
    Ts = Ts * scaling;

    %These plots are useful just for comparing that the previous velocities
    %multiplied by the scaling factor will be equal to the new velocities
    %obtained with the new timing law.
    figure;
    subplot(2,1,1);
    plot(t_scaled, v_t * scaling); % Scaled linear velocity
    xlabel('Time (s)');
    ylabel('Linear Velocity (m/s)');
    title('Scaled Linear Velocity');

    %Plot angular velocity
    subplot(2,1,2);
    plot(t_scaled, w_t * scaling); % Scaled angular velocity
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Scaled Angular Velocity');
    hold off;

    % Repeat the computation for computing the path with the scaled time
    trajDuration = max(t_scaled);
    a0 = 0;
    a1 = 0;
    a2 = -3 * a0 / trajDuration^2 - 2 * a1 / trajDuration + 3 / trajDuration^2 - 0 / trajDuration;
    a3 = 2 * a0 / trajDuration^3 + a1 / trajDuration^2 - 2 / trajDuration^3 + 0 / trajDuration^2;

    for i=1:length(t_scaled)
       [s_scaled(i), ds_scaled(i), dds_scaled(i)] = cubic_polinomial(a0, a1, a2, a3, t_scaled(i));
    end
    x_s_scaled = s_scaled.^3 * x_f - (s_scaled - 1).^3 * x_i + alpha_x * s_scaled.^2 .* (s_scaled - 1) + beta_x * s_scaled .* (s_scaled - 1).^2;
    y_s_scaled = s_scaled.^3 * y_f - (s_scaled - 1).^3 * y_i + alpha_y * s_scaled.^2 .* (s_scaled - 1) + beta_y * s_scaled .* (s_scaled - 1).^2;
    x_s=x_s_scaled; 
    y_s=y_s_scaled;
    % Renaming the trajectories will be useful for the 4th exercise, in
    % such way that trajectory without scaling and with scaling will have
    % the same name
    % Using exactly the same name would have been a solution to make code
    % smaller, but it's for a better clarification of actions.
    
    x_p_scaled = 3 * x_f * s_scaled.^2 - 3 * x_i * (s_scaled - 1).^2 + alpha_x * s_scaled.^2 + 2 * alpha_x * (s_scaled - 1) .* s_scaled + beta_x * (s_scaled - 1).^2 + 2 * beta_x * s_scaled .* (s_scaled - 1);
    y_p_scaled = 3 * y_f * s_scaled.^2 - 3 * y_i * (s_scaled - 1).^2 + alpha_y * s_scaled.^2 + 2 * alpha_y * (s_scaled - 1) .* s_scaled + beta_y * (s_scaled - 1).^2 + 2 * beta_y * s_scaled .* (s_scaled - 1);
    x_p=x_p_scaled;
    y_p=y_p_scaled;

    x_pp_scaled = 6 * x_f * s_scaled - 3 * x_i * (s_scaled - 1) + 2 * alpha_x * (s_scaled - 1) + 4 * alpha_x * s_scaled + 4 * beta_x * (s_scaled - 1) + 2 * beta_x * s_scaled;
    y_pp_scaled = 6 * y_f * s_scaled - 3 * y_i * (s_scaled - 1) + 2 * alpha_y * (s_scaled - 1) + 4 * alpha_y * s_scaled + 4 * beta_y * (s_scaled - 1) + 2 * beta_y * s_scaled;
    x_pp=x_pp_scaled;
    y_pp=y_pp_scaled;

    v_s_scaled = sqrt(x_p_scaled.^2 + y_p_scaled.^2);
    w_s_scaled = (y_pp_scaled .* x_p_scaled - x_pp_scaled .* y_p_scaled) ./ (x_p_scaled.^2 + y_p_scaled.^2);
    v_s=v_s_scaled;
    w_s=w_s_scaled;

    v_t_scaled = v_s_scaled .* ds_scaled;
    w_t_scaled = w_s_scaled .* ds_scaled;
    v_t = v_t_scaled;
    w_t = w_t_scaled;

    %The equalities after the computation of the scaled quantities
    theta = atan2(y_p_scaled,x_p_scaled);
    figure
    plot(t_scaled,theta);
    xlabel('Times (s)');
    ylabel('Theta');
    title('Robot angular position');

    % Plot della traiettoria scalata
    figure;
    plot(x_s_scaled, y_s_scaled, 'LineWidth', 2);   % The path computed with scaled trajectories will remain the same as before
    hold on;
    plot(qi(1), qi(2), 'ro', 'MarkerSize', 10); % Plot initial configuration
    plot(qf(1), qf(2), 'go', 'MarkerSize', 10); % Plot final configuration
    xlabel('x (m)');
    ylabel('y (m)');
    title('Robot Path (Scaled Time)');
    legend('Path', 'qi', 'qf');
    axis equal;
    hold off;

    % Plot delle velocità scalate
    figure;
    subplot(2,1,1);
    plot(t_scaled, v_s);
    xlabel('Time (s)');
    ylabel('Linear Velocity (m/s)');
    title('Linear Velocity v(s) (Scaled Time)');
    subplot(2,1,2);
    plot(t_scaled, w_s);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity w(s) (Scaled Time)');

    figure
    subplot(2,1,1)
    plot(t_scaled, v_t_scaled);
    xlabel('Time (s)');
    ylabel('Linear Velocity (m/s)');
    title('Linear Velocity v(t) (Scaled Time)');
    subplot(2,1,2);
    plot(t_scaled, w_t_scaled);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity  w(t) (Scaled Time)');
else
    % If velocity constraints are satisfied, proceed without scaling time
    disp('Velocity constraints are satisfied. No need to scale time.');
end

%% Exercise n.4

b=0.2; %Let%s choose a point B at a distance |b| from the center of the wheel along the sagittal axis.
       %b>0 means that the point is ahead the unicycle

t_sim=max(t);

%Choose k1,k2>0
k1=5;
k2=5;

x_t=timeseries(x_s,t);
y_t=timeseries(y_s,t);
v_t=timeseries(v_t,t);
w_t=timeseries(w_t,t);
theta_t=timeseries(theta,t);

syms th;    %theta
syms bb;    %b
T=[cos(th), -bb*sin(th);sin(th),bb*cos(th)];
T_inv=inv(T);       % The computation of T using syms variables helps to 
                    % compute the inverse of T to be put in the
                    % Simulink model
simplify(T_inv);

%% Function
% Cubic polynomial function
function [s, ds, dds] = cubic_polinomial(a0, a1, a2, a3, t)
    s = a3 * t.^3 + a2 * t.^2 + a1 * t + a0;
    ds = 3 * a3 * t.^2 + 2 * a2 * t + a1;
    dds = 6 * a3 * t + 2 * a2;
end