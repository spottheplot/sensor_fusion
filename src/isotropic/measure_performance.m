function [RMSE_EKF, CRLB_EKF, RMSE_UKF, CRLB_UKF, RMSE_PF, CRLB_PF] = measure_performance( plotting, nIterations)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
% times, it will calculate the RMSE and plot the results in svg format.
% This function is temporary so no effort in proper and efficient coding is done.
% Example of use [RMSE_EKF, CRLB_EKF, RMSE_UKF, CRLB_UKF, RMSE_PF, CRLB_PF] = measure_performance(false, 10)
% TODO: This function should be run after tuning individual Q and R for
% each filter.
% TODO2: Add possibility of modifying P_init


%% Initialisation of variables for each filter
% Matrix for storing each state (2) at each increment of time (1801) for
% every iteraion (nIterations)
x_state_ekf = zeros(2, 1801, nIterations);
x_state_ukf = zeros(2, 1801, nIterations);
x_state_pf = zeros(2, 1801, nIterations);
% RMSE for each filter and for each t of all the iterations
RMSE_EKF = zeros(nIterations, 1801);
CRLB_EKF = zeros(nIterations, 1801);
RMSE_UKF = zeros(nIterations, 1801);
CRLB_UKF = zeros(nIterations, 1801);
RMSE_PF = zeros(nIterations, 1801);
CRLB_PF = zeros(nIterations, 1801);

Q = 0.1;
R = 0.05;

for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

for i=1:nIterations
    [x_state_ekf(:, :, i), x_t_vec, P_cov(:,:,:,i)] = Main_isotropic_EKF(plotting, Q, R, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i));
    % For each t calculate the distance between the estimate and the real location
    for k=1 : size(x_state_ekf, 2)
        % Calculate the norm of all the states 
        RMSE_EKF(i, k) = norm((x_state_ekf(:,k, i)- x_t_vec'));
        CRLB_EKF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
% Calculate the average error along all iterations for each t in the process
RMSE_EKF = mean(RMSE_EKF);
CRLB_EKF = mean(CRLB_EKF);

for i=1:nIterations
    [x_state_ukf(:, :, i), x_t_vec, P_cov(:,:,:,i)] = Main_isotropic_UKF(plotting, Q, R, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i));
    for k=1 : size(x_state_ukf, 2)
        RMSE_UKF(i, k) = norm((x_state_ukf(:,k, i)- x_t_vec'));
        CRLB_UKF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_UKF = mean(RMSE_UKF);
CRLB_UKF = mean(CRLB_UKF);

for i=1:nIterations
    [x_state_pf(:, :, i), x_t_vec, P_cov(:,:,:,i)] = Main_isotropic_PF(plotting, 10*Q, 10*R, x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i));
    for k=1 : size(x_state_pf, 2)
        RMSE_PF(i, k) = norm((x_state_pf(:,k, i)- x_t_vec'));
        CRLB_PF(i, k) = norm(diag(sqrt(diag(P_cov(:,:,k,i))))); % Eliminate the cross terms of the inverse of the covariance matrix
    end;
end
RMSE_PF = mean(RMSE_PF);
CRLB_PF = mean(CRLB_PF);

%% Individual plots for each filter -- Filter vs Cramer Rao Lower Bound
mkdir('performance_comparison')
% EKF
CRLB_EKF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_EKF)
hold on
plot(CRLB_EKF, 'g')
legend('RMSE', 'CRLB')
title('EKF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./performance_comparison/EKF','-dsvg')
% EKF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track EKF RMSE (m)')
print('./performance_comparison/EKF_final','-dsvg')
% UKF
CRLB_UKF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_UKF)
hold on
plot(CRLB_UKF, 'g')
legend('RMSE', 'CRLB')
title('UKF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./performance_comparison/UKF','-dsvg')
% UKF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track UKF RMSE (m)')
print('./performance_comparison/UKF_final','-dsvg')

% PF
CRLB_PF(1,1:25) = NaN;
figure('Visible','off')
plot(RMSE_PF)
hold on
plot(CRLB_PF, 'g')
legend('RMSE', 'CRLB')
title('PF RMSE (m)')
xlabel('Time step')
ylabel('RMSE (m)')
print('./performance_comparison/PF','-dsvg')
% PF-Final
xlim([1000, 1800])
ylim('auto')
title('Final track PF RMSE (m)')
print('./performance_comparison/PF_final','-dsvg')

%% Comparison plots between filters -- EKF vs UKF vs PF
clear plot;
figure('Visible','off');
plot(RMSE_EKF, 'r');
hold on;
plot(RMSE_UKF, 'b');
plot(RMSE_PF , 'g');
legend ('EKF','UKF','PF')
xlabel('Time step')
ylabel('RMSE (m)')
print('./performance_comparison/All','-dsvg')
% Final part
xlim([1000, 1800])
ylim('auto')
title('Final track (m)')
print('./performance_comparison/All_final','-dsvg')
end

