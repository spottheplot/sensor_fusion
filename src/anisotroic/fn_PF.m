function [X, P_f] = fn_PF(uav_init_pos, uav_actual_pos, h0, alpha_meas, F, G, Q, R)
persistent particles;

nsamples =500;

if isempty(particles)
   particles(1,:)=4000*rand(1,nsamples) + 4000; 
   particles(2,:)=4000*rand(1,nsamples) + 4000;
end

% Evaluate weights:
pred_state = zeros(size(particles, 1), nsamples);
pred_meas = zeros(size(alpha_meas, 2), nsamples);
for i = 1 : nsamples
pred_state(:,i)   = fk(particles(:,i));      % predicted measurements from samples
pred_meas(:,i) = hk(pred_state(:,i));
end

likelihood=zeros(1,nsamples);
for i = 1 : size(pred_meas, 2)
    likelihood(i) = 1/sqrt((2*pi*det(R)))*exp(-0.5*((alpha_meas - pred_meas(:, i))'/R*(alpha_meas - pred_meas(:, i))));    % evaluate likelihoods 
end
weight=cumsum(likelihood/sum(likelihood));      % normalise weights & form cumulative distribution    
%  re-sampling procedure (systematic)    
addit=1/nsamples; 
stt=addit*rand(1);     
selection_points=[ stt : addit : stt+(nsamples-1)*addit ];

j=1;       %set up comb 
x_post=zeros(size(pred_state));
for i=1:nsamples  
    while selection_points(i) >= weight(j);        
        j=j+1;      
    end;     
    x_post(:,i)=pred_state(:,j);              
end;

% figure(2);
% % plot3(pred_state(1,:), pred_state(2,:), likelihood,'+');
% hold on;
% plot(pred_state(1,:),pred_state(2,:),'.b','markersize',30);
% plot(x_post(1,:),x_post(2,:), '.m','markersize',10);
% pause;
% close;
% figure(1);

particles = x_post;

X = mean(particles');

P_f = diag(var(particles, 0,2));

%         samp_store(:,:,k) = x_post;        % output:  store posterior samples (for analysis only)

       
%% h(X): Nonlinear measurement eq
function h=hk(x_particle) %TODO
    alpha_clean = norm([x_particle; 0] - [uav_init_pos, h0]')^2 / norm([x_particle; 0] - [uav_actual_pos, h0]')^2;
    h = alpha_clean; %+ randn(size(R,2), size(alpha_clean,1)) * diag(sqrt(R));
end
%% partial_h/partial_X: Jacobian of measurement eq
function f=fk(x_particle) %TODO
    f = F * x_particle + G * sqrt(Q) * randn(2, 1);
end

end