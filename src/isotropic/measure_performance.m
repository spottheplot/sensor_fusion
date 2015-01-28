function [RMSE_EKF, RMSE_UKF, RMSE_PF] = measure_performance( plot, nIterations)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
%times and will calculate the RMSE.
%   This function is temporary so no effort in proper and efficient coding is done.

x_state_ekf = zeros(2, 1801, nIterations);
x_state_ukf = zeros(2, 1801, nIterations);
x_state_pf = zeros(2, 1801, nIterations);
RMSE_EKF = zeros(nIterations, 1801);
RMSE_UKF = zeros(nIterations, 1801);
RMSE_PF = zeros(nIterations, 1801);

Q = 0.1;
R = 0.05;

for iteration=1:nIterations
  x_jammer(:,:,iteration) = place_jammer();
  [x_uav(:,:,iteration), psi_uav(:,:,iteration)] = place_uav();
end

for iteration=1:nIterations
    [x_state_ekf(:, :, iteration), x_t_vec] = Main_isotropic_EKF(false, Q, R, x_jammer(:,:,iteration), x_uav(:,:,iteration), psi_uav(:,:,iteration));
    for l=1 : size(x_state_ekf, 2)
        RMSE_EKF(iteration, l) = norm((x_state_ekf(:,l, iteration)- x_t_vec'));
    end;
end
% plot(1:1:size(RMSE_EKF, 2), mean(RMSE_EKF))
RMSE_EKF = mean(RMSE_EKF);

for iteration=1:nIterations
    [x_state_ukf(:, :, iteration), x_t_vec] = Main_isotropic_UKF(false, Q, R, x_jammer(:,:,iteration), x_uav(:,:,iteration), psi_uav(:,:,iteration));
    for l=1 : size(x_state_ukf, 2)
        RMSE_UKF(iteration, l) = norm((x_state_ukf(:,l, iteration)- x_t_vec'));
    end;
end
% plot(1:1:size(RMSE_UFK, 2), mean(RMSE_UFK))
RMSE_UKF = mean(RMSE_UKF);

% for iteration=1:nIterations
%     [x_state_pf(:, :, iteration), x_t_vec] = Main_isotropic_PF(false, Q, R, x_jammer(:,:,iteration), x_uav(:,:,iteration), psi_uav(:,:,iteration));
%     for l=1 : size(x_state_pf, 2)
%         RMSE_PF(iteration, l) = norm((x_state_pf(:,l, iteration)- x_t_vec'));
%     end;
% end
% % plot(1:1:size(RMSE_PF, 2), mean(RMSE_PF))
% RMSE_PF = mean(RMSE_PF);

%%
% for iteration=1:nIterations
%     [x_state_ekf(:, :, iteration), x_t_vec] = Main_isotropic_EKF(false, x_jammer, x_uav);
%     for l=1 : size(x_state_ekf, 2)
%         RMSE_EKF(iteration, l) = norm((x_state_ekf(:,l, iteration)- x_t_vec'));
%     end;
% end
% % plot(1:1:size(RMSE_EKF, 2), mean(RMSE_EKF))
% RMSE_EKF = mean(RMSE_EKF);
% 
% hold on
% 
% for iteration=1:nIterations
%     [x_state_ukf(:, :, iteration), x_t_vec] = Main_isotropic_UKF(false, x_jammer);
%     for l=1 : size(x_state_ukf, 2)
%         RMSE_UFK(iteration, l) = norm((x_state_ukf(:,l, iteration)- x_t_vec'));
%     end;
% end
% % plot(1:1:size(RMSE_UFK, 2), mean(RMSE_UFK))
% RMSE_UKF = mean(RMSE_UKF);
% 
% for iteration=1:nIterations
%     [x_state_pf(:, :, iteration), x_t_vec] = Main_isotropic_EKF(false, x_jammer);
%     for l=1 : size(x_state_pf, 2)
%         RMSE_PF(iteration, l) = norm((x_state_pf(:,l, iteration)- x_t_vec'));
%     end;
% end
% % plot(1:1:size(RMSE_PF, 2), mean(RMSE_PF))
% RMSE_PF = mean(RMSE_PF);

end

