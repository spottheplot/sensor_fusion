function [meanRMSEvalues, meanRMSE_finalValues] = tune_filter( plotting, nIterations, filter_function)
%PERFORMANCE_CHECKER This function will run the three filters nIterations
%times and will calculate the RMSE.
%  This function is temporary so no effort in proper and efficient coding is done.
%  INPUT
% plotting - Boolean to determine if we want to plot the proccess while simulating
% nIterations - Number of iterations to run
% filter_function - A function handler pointing to the filter rhat should
% be used.
% EXAMPLE OF USE
% [meanRMSEvalues, meanRMSE_finalValues] = tune_filter( false, 10, @Main_isotropic_EKF)
% TODO: Allow tuning of  P_init (Initial covariance matrix)

%% Initialise variables
x_state = zeros(2, 1801, nIterations);
RMSE = zeros(nIterations, 1801);
% Values over we want to iterate
Q = [0.001, 0.01, 0.1, 1, 10, 100];
R = [0.001, 0.01, 0.1, 1, 10, 100];

meanRMSEvalues = zeros(size(Q,2), size(R,2));
meanRMSE_finalValues = zeros(size(Q,2), size(R,2));

% Declaring function to get
%% Start iterations
for i=1:nIterations
  x_jammer(:,:,i) = place_jammer();
  [x_uav(:,:,i), psi_uav(:,:,i)] = place_uav();
end

for m=1:size(Q,2)
    for k=1:size(R,2)
        for i=1:nIterations
            [x_state(:, :, i), x_t_vec] = filter_function(plotting, Q(m), R(k), x_jammer(:,:,i), x_uav(:,:,i), psi_uav(:,:,i));
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
surf(Q,R, meanRMSEvalues)
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10)
xlabel('Q')
ylabel('R')
zlabel('RMSE')
shading interp
print(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations)),'-dsvg')
% Plotting only RMSE of final values to evaluate convergence quality
h = figure;
surf(Q,R, meanRMSE_finalValues)
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10)
xlabel('Q')
ylabel('R')
zlabel('RMSE')
shading interp
print(h, strcat('./tuning_results/', func2str(filter_function), int2str(nIterations), '_FINAL'),'-dsvg')
end

