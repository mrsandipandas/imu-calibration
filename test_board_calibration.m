%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU-IMU calibration code with sample logs collected on the test
% board. The logs are provided for quick evaluation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
set(0,'defaultTextInterpreter','latex'); %Set the default for rendering
legend_font_size = 20;
label_font_size = 24;
tick_font_size = 18;
line_width = 2;
title_size = 30;
axis_desc_size = 25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test different logs collected with the IMU board
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Global variables for the generated solution
global x_star  err_iter

% Change the log files name to test other configurations
file_a = 'logs/a_90_2.mat'; 
file_b = 'logs/b_90_2.mat';

% IMU-B is origin. IMU-A is rotated and its parameters are estimated wrt IMU-B. 
% The GT is manually measured with a tape.
t_a = [-0.19 0.197 0.00]; % x y z position
R_a = eul2rotm([0.00 0.00 -pi/2],'XYZ'); % roll, pitch, yaw
T_a = [R_a t_a'; zeros(1,3) 1];

% Configuration parameters
use_solution = 2; %Which solution we use for visualization  1: Closed form LS, 2: BVLS
fs = 100; %IMU frequency
n_IMU = 2;
bound_config = 0.05; % (in metres)
smoothen_factor = 10; % How many samples to add and then average for smoothening out signal
G = 9.83;
window_size = 50; %Window size for information matrix calculation
obs_thershold = 9e-14; % Thershold for observability, this needs to be tuned

% Transformation matrices, based on the test board setup
t_b = [0 0 0]; % x y z position
R_b = eul2rotm([0.00 0.00 0.00],'XYZ'); % roll, pitch, yaw
T_b = [R_b t_b'; zeros(1,3) 1];

T_star = T_b \ T_a;
x_star = T_star(1:3,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data and concatenate
% data format for storing [acc(imu1/imu2/imu3/...) omega (imu1/imu2/imu3/...)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfile(file_a) && isfile(file_b)    
    N = n_IMU;
    
    % get the size of data, the sizes are different in each IMU, get min
    M = inf;
    for i = 1: N
        if (i<N)
            data = load(file_a);
            imu_field = strcat('IMU_','a');
        else
            data = load(file_b);
            imu_field = strcat('IMU_','b');
        end
        M = min(M, length(data.(imu_field).Acc_X));
    end

    % down sample the data
    M = floor(M/smoothen_factor);

    % initializations
    acc = zeros(N,M,3);
    gyro = zeros(N,M,3);
    gyro_acc = zeros(N,M,3); %Note the first data will be junk, Please omit it
    A = zeros(N,3*(M-1),3);
    B = zeros(3*(M-1),1);
    %----------------------------------------------------------------------
    for i = 1:N % loop over all IMUs
        if (i<N)
            data = load(file_a);
            imu_field = strcat('IMU_','a');
        else
            data = load(file_b);
            imu_field = strcat('IMU_','b');
        end
        
        acc_tmp = [0 0 0];
        gyro_tmp = [0 0 0];
        prev_gyro_tmp = [0 0 0];
        prev_time_tmp = 0.0;
        time_tmp = 0.0;
        A_i = [];
        
        for k = 1:M*smoothen_factor % loop over all data
            acc_tmp = acc_tmp + [data.(imu_field).Acc_X(k) data.(imu_field).Acc_Y(k) data.(imu_field).Acc_Z(k)] - [0 0 G];
            gyro_tmp = gyro_tmp + [data.(imu_field).Gyr_X(k) data.(imu_field).Gyr_Y(k) data.(imu_field).Gyr_Z(k)];
            time_tmp = double(data.(imu_field).SampleTimeFine(k))/10000;
                
            if (mod(k,smoothen_factor) == 0)
                index = k/smoothen_factor;
                
                % read IMU data
                acc(i,index,:) = acc_tmp/smoothen_factor;

                gyro_tmp = gyro_tmp/smoothen_factor;
                
                gyro(i,index,:) = gyro_tmp;
                
                gyro_acc_tmp = (gyro_tmp - prev_gyro_tmp)/(time_tmp - prev_time_tmp);
                
                gyro_acc(i,index,:) = gyro_acc_tmp;
                
                if (index > 1) %Skip the first one as the angular acc is not valid
                    A_i = [A_i; getSkewMatrix(gyro_tmp)*getSkewMatrix(gyro_tmp) + getSkewMatrix(gyro_acc_tmp)];
                end

                prev_time_tmp = time_tmp;
                prev_gyro_tmp = gyro_tmp;
                acc_tmp = [0 0 0];
                gyro_tmp = [0 0 0];
                time_tmp = 0;
            end
        end
        A(i,:,:) = A_i;
    end

% To access data in Mx3 format
acc_a = reshape(acc(1,:,:),[],3);
acc_b = reshape(acc(2,:,:),[],3);   
gyro_a = reshape(gyro(1,:,:),[],3);
gyro_acc_a = reshape(gyro_acc(1,:,:),[],3);
gyro_b = reshape(gyro(2,:,:),[],3);
gyro_acc_b = reshape(gyro_acc(2,:,:),[],3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibration estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 1/fs;

% Do noise compensation according to the paper or other way if required. We
% left it out for simplicity and focused mainly on the calibration part.

% Estimating R with Umeyama alignment
R = computeRotation(gyro_a', gyro_b');

% Estimating t
B_i = [];
for i = 1:M
if (i>1) %Skip the first one as the angular acc is not valid
    B_i = [B_i; R'*acc_b(i,:)' - acc_a(i,:)'];
end
end
B = B_i;

A = reshape(A(1,:,:),[],3);

t_GT = x_star'
t_closed_form = ((A'*A)\A'*B)'


x0 = t_a'; % Give the positional values
bound = bound_config*[1 1 1]';

err_iter = [];

% Algorithm: 'active-set', 'interior-point', 'sqp', 'trust-region-reflective', or 'sqp-legacy'
options = optimoptions('fmincon','Algorithm','sqp', 'Display','iter', 'OutputFcn',@outfun);
options.OptimalityTolerance = 1e-10;
options.ConstraintTolerance = 1e-10;

fun = @(x)getObjectiveFunction(A, x, B);
[x,fval] = fmincon(fun,x0,[],[],[],[],(x0 - bound),(x0 + bound),[],options);

t_bvls = x'

R_GT = rotm2eul(T_star(1:3,1:3), 'XYZ')*180/pi
R_estimated = rotm2eul(R, 'XYZ')*180/pi

hold on; grid on;
plot(err_iter, '--.');
title('optimization error norm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch use_solution
    case 1 %Closed form LS
        t_hat = t_closed_form';
        R_hat = R;
    case 2 %BVLS
        t_hat = x;
        R_hat = R;
    otherwise %No calibration
        t_hat = [0 0 0]';
        R_hat = eye(3);
end

omega_i = [];
for i = 1:M
    omega_i = [omega_i; (R_hat*(getSkewMatrix(gyro_a(i,:))*getSkewMatrix(gyro_a(i,:)) + getSkewMatrix(gyro_acc_a(i,:)))*x)'];
end

plotResults('Before Calibration', acc_a, acc_b, gyro_a, gyro_b, smoothen_factor, legend_font_size, label_font_size, tick_font_size, line_width, title_size, axis_desc_size)

acc_a = (R_hat*acc_a')' + omega_i;
gyro_a = (R_hat*gyro_a')';

plotResults('After Calibration', acc_a, acc_b, gyro_a, gyro_b, smoothen_factor, legend_font_size, label_font_size, tick_font_size, line_width, title_size, axis_desc_size)

% Observability analysis
isValid = zeros(M,1);
for i=1:window_size:ceil(M/window_size)*window_size
    if (i+window_size-1 > M)
        index = M;
    else
        index = i+window_size-1;
    end

    gyro_a_batch = gyro_a(i:index,:);
    gyro_b_batch = gyro_b(i:index,:);


    [U,S,V] = svd(gyro_a_batch*inv(cov((R*gyro_a_batch'-gyro_b_batch')'))*gyro_a_batch');
    
    if (min(diag(S)) > obs_thershold)
        isValid(i:index) = 1;
    end
end
plotObservabilityAnalysis (gyro_a, isValid, window_size, legend_font_size, label_font_size, tick_font_size, line_width, title_size, axis_desc_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = computeRotation(A, B) %Umeyama Alignment
    [~,n]=size(A);

    %Mean Center Data
    Ac = mean(A')';
    Bc = mean(B')';
    A = A-repmat(Ac,1,n);
    B = B-repmat(Bc,1,n);

    %Calculate Optimal Rotation
    [u,~,v]=svd(A*B');
    R = v*diag([1 1 det(v*u')])*u';
    %t = Bc - R*Ac;
end

function plotResults(plot_description, acc_a, acc_b, gyro_a, gyro_b, avg_factor, legend_font_size, label_font_size, tick_font_size, line_width, title_size, axis_desc_size)
    %%%%%%%%%%%%%%%%%%%%%% Acc data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('name', strcat('IMU Data (', plot_description , ')'));
    subplot(3,2,1);
    plot(acc_a(:,1), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    hold on
    plot(acc_b(:,1), 'Linewidth', line_width)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('m/s${}^2$', 'FontSize', label_font_size)
    grid on
    title('$a_x$', 'interpreter','latex', 'FontSize', title_size)

    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)

    subplot(3,2,3);
    plot(acc_a(:,2), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    hold on
    plot(acc_b(:,2), 'Linewidth', line_width)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('m/s${}^2$', 'FontSize', label_font_size)
    grid on
    title('$a_y$', 'interpreter','latex', 'FontSize', title_size )
    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)

    subplot(3,2,5);
    plot(acc_a(:,3), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    hold on
    plot(acc_b(:,3), 'Linewidth', line_width)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('m/s${}^2$', 'FontSize', label_font_size)
    grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%title('$a_z$', 'interpreter','latex', 'FontSize', title_size )
    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)

    %%%%%%%%%%%%%%%%%%%%%% Gyro data %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(3,2,2);
    plot(gyro_a(:,1), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    hold on
    plot(gyro_b(:,1), 'Linewidth', line_width)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on
    title('$\omega_x$', 'interpreter','latex', 'FontSize', title_size )
    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)

    subplot(3,2,4);
    plot(gyro_a(:,2), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    hold on
    plot(gyro_b(:,2), 'Linewidth', line_width)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on
    title('$\omega_y$', 'interpreter','latex', 'FontSize', title_size )
    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)

    subplot(3,2,6);
    plot(gyro_a(:,3), 'Linewidth', line_width)
    ax=gca;
    ax.FontSize = tick_font_size;
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/avg_factor)
    hold on
    plot(gyro_b(:,3), 'Linewidth', line_width)
    xlabel('Time (seconds)', 'FontSize', label_font_size)
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on
    title('$\omega_z$', 'interpreter','latex', 'FontSize', title_size )
    legend('$\mathtt{IMU}-\mathtt{A}$', '$\mathtt{IMU}-\mathtt{B}$', 'interpreter','latex', 'FontSize',legend_font_size)
end

% return skew-symmetric matrix given a vector a (3x1)
function m = getSkewMatrix(a)
m = [
    0    , -a(3), a(2)
    a(3) , 0    , -a(1)
    -a(2), a(1) , 0 ]; %3x3
end

% objective function for optimizatoin
function [c,ceq] = getObjectiveFunction(A, x, B)
    ceq = [];
    c = (norm(A*x - B))^2;
end

% output function to get iteration variables
function stop = outfun(x,optimValues,state)
global x_star err_iter
stop = false;
switch state
    case 'init'
    case 'iter'
        err = norm(x-x_star)/norm(x_star);
        err_iter = [err_iter err];
    case 'done'
    otherwise
end
end

function drawPatch(isValid, window_size)
    for i=1:window_size:ceil(length(isValid)/window_size)*window_size
        if (i+window_size-1 > length(isValid))
            index = length(isValid);
        else
            index = i+window_size-1;
        end
        if (isValid(i) == 1 && isValid(index))
            patch([i index index i], [min(ylim) min(ylim) max(ylim) max(ylim)], 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);
        end
    end
end

function plotObservabilityAnalysis (gyro_data, isValid, window_size, legend_font_size, label_font_size, tick_font_size, line_width, title_size, axis_desc_size)
%%%%%%%%%%%%%%%%%%%%%% Plot all Gyro data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('name','Observability analysis');
    subplot(3,1,1);
   
    plot(gyro_data(:,1))
    hold on
    %plot(gyro_FRT_original(:,1))
    drawPatch(isValid, window_size);
    xlabel('Time (seconds)', 'FontSize', label_font_size);
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on
    box off
    title('$\omega_x$', 'fontsize',30,'interpreter','latex')
    %legend('FLT')
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/10)

    a = get(gca,'XTickLabel');
    
    set(gca,'XTickLabel',a,'fontsize',24)
    hold off

    subplot(3,1,2);
    plot(gyro_data(:,2))
    hold on
    %plot(gyro_FRT_original(:,2))
    drawPatch(isValid, window_size);
    xlabel('Time (seconds)', 'FontSize', label_font_size);
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on; box off;
    title('$\omega_y$', 'fontsize',30,'interpreter','latex')
    %legend('FLT')
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/10)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',24)
    hold off

    subplot(3,1,3);
    plot(gyro_data(:,3))
    hold on
    %plot(gyro_FRT_original(:,3))
    drawPatch(isValid, window_size);
    xlabel('Time (seconds)', 'FontSize', label_font_size);
    ylabel('rad/s', 'FontSize', label_font_size)
    grid on
    box off
    title('$\omega_z$', 'fontsize',30,'interpreter','latex')
    %legend('FLT')
    xt = get(gca,'xtick');
    set(gca,'XTick',xt, 'xticklabel',xt/10)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',24)
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

