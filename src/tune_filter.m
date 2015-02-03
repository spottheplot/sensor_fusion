function [meanRMSEvalues, meanRMSE_finalValues, Q, R, P_init] = tune_filter( plotting, nIterations, filter_function)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
%times and will calculate the RMSE.
%  This function is temporary so no effort in proper and efficient coding is done.
%  INPUT
% plotting - Boolean to determine if we want to plot the proccess while simulating
% nIterations - Number of iterations to run
% filter_function - A function handler pointing to the filter rhat should
% be used.
% OUTPUT 
% meanRMSEvalues Is the mean of the RMSE of all the iterations. Used to measure general performance
% meanRMSE_finalValues is the mean of the RMSE of the last iterations. Used to measure the quality of the convergence.
% EXAMPLE OF USE
% [meanRMSEvalues, meanRMSE_finalValues] = tune_filter( false, 10, @Main_isotropic_EKF)
% TODO: Allow tuning of  P_init (Initial covariance matrix)

%% Initialise variables
x_state = zeros(2, 1801, nIterations);
RMSE = zeros(nIterations, 1801);
% Values over we want to iterate
Q = [0.001, 0.01, 0.1, 1, 10];
R = [0.01, 0.1, 1, 10];
P_init = 400;

meanRMSEvalues = zeros(size(Q,2), size(R,2));
meanRMSE_finalValues = zeros(size(Q,2), size(R,2));

% Global variables needed for place_jammer() and place_uav()
global x_bnd y_bnd d2r
x_bnd=12*10^3;                                                      %   x area boundary [m]
y_bnd=12*10^3;                                                      %   y area boundary [m]
d2r=pi/180;                                                         %   Value in rad = Value in deg * d2r
% Jammer and UAV variables
x_jammer = zeros(1,2,nIterations);
x_uav = zeros(1,2,nIterations);
psi_uav = zeros(1,1,nIterations);

%% Start iterations
for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

for m=1:size(Q,2)
    for k=1:size(R,2)
        for i=1:nIterations
            [x_state(:, :, i), x_t_vec] = filter_function(plotting, Q(m), R(k), x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i),P_init);
            for l=1 : size(x_state, 2)
                RMSE(i, l) = norm((x_state(:,l, i)- x_t_vec'));
            end
        end
        meanRMSEvalues(m,k) = mean(mean(RMSE));
        meanRMSE_finalValues(m,k) = mean(mean(RMSE(:, 1500:end)));
        RMSE = zeros(nIterations, 1801);
    end
end

mkdir('tuning_results')
h = figure;
surf(Q,R, meanRMSEvalues')
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
caxis([0 2000])
colormap jet
xlabel('Q')
ylabel('R')
zlabel('RMSE (m)')
title(sprintf('RMSE vs Q & R.\n P_init = %d', P_init))
shading interp
print(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations)),'-dsvg')
savefig(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations)))
% Plotting only RMSE of final values to evaluate convergence quality
h = figure;
surf(Q,R, meanRMSE_finalValues')
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZLim', [0 2000], 'FontSize', 10)
caxis([0 2000])
colormap jet
shading interp
xlabel('Q')
ylabel('R')
zlabel('RMSE (m)')
title(sprintf('Final track RMSE vs Q & R.\n Pinit: %d^2 Min: %d m', P_init, round(min(min(meanRMSE_finalValues)))))
print(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations), '_FINAL'),'-dsvg')
savefig(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations), '_FINAL'))
end

