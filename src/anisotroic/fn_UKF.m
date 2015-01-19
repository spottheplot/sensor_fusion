function [X, P_f, K] = fn_UKF(uav_init_pos, uav_actual_pos, h0, alpha_meas, X_init, P_init, F, G, Q, R)
%% Unscented Kalman Filter
%% =========================
%% Prediction
%% =========================
%% Augmenting state and covariance with mean and covariance of process noise
% X_init = [X_init; 0; 0];
% P_init = [P_init zeros(size(P_init,2),size(Q,1)); zeros(size(Q,1), size(P_init,2)) Q];
%% Parameters
L=numel(X_init); 
alpha=1e-4;                                 %default, tunable
ki=1;                                       %default, tunable
beta=0.5;                                   %default, tunable
lambda=alpha^2*(L+ki)-L;                    %scaling factor
c=L+lambda;                                 %scaling factor
%% Create set of sigma points
X=sigmas(X_init, P_init, c);
%% Create specifig weights
Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
%% Project points through dynamic function
[x1,X1,P1,X2]=ut(@fk,X,Wm,Wc,L,Q)  ;        %unscented transformation of process
%% =========================
%% Correction
%% =========================
%% Augmenting with mean and covariance of measurement noise
% x1 = [x1; 0; 0];
% P1 = [P1 zeros(size(P1,2),size(R,1)); zeros(size(R,1), size(P1,2)) R];
% %% Parameters
% L=numel(x1); 
% alpha=1e-3;                                 %default, tunable
% ki=0;                                       %default, tunable
% beta=2;                                     %default, tunable
% lambda=alpha^2*(L+ki)-L;                    %scaling factor
% c=L+lambda;                                 %scaling factor
%% Create set of sigma points
% X1=sigmas(x1,P1,c);                         %sigma points around x1
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
%% Create specifig weights
% Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
% Wc=Wm;
% Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
%% Project points through observation function
[z1,~,P2,Z2]=ut(@hk,X1,Wm,Wc,1,R);          %unscented transformation of measurments
P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
%% Calculate K
K=P12/(P2);
%% Update state and covariance
X=x1+K*(alpha_meas-z1);                              %state update
P_f=P1-K*P12';                                %covariance update
%% ===============================================

%% h(X): Nonlinear measurement eq
function h=hk(X_predicted)
h=norm([X_predicted; 0] - [uav_init_pos, h0]')^2 / norm([X_predicted; 0] - [uav_actual_pos, h0]')^2;
end

%% h(X): Dynamic equation
function f=fk(X_predicted)
f = X_predicted;
end

function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2); % Number of sigma points
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L
    Y(:,k)=f(X(:,k));       
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1' + R;
end

end

function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

% A = c*chol(P)';
A = sqrt(c)*chol(P)';
Y = x(:,ones(1,numel(x))); % Creates L points;
X = [x Y+A Y-A]; % Creates 2L + 1 points
end