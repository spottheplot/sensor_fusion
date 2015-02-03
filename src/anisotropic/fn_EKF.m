function [X, P_f, K] = fn_EKF(uav_init_pos, uav_actual_pos, h0, alpha_meas, X_k, P_k, F, G, Q, R)
%% KF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF);

%% Extended Kalman Filter
%% =========================
%% Prediction
%% =========================

%% Equation 1
X_s = F * X_k;
%% Equation 2
P_s = F * P_k * F' + G * Q * G'; % Error covariance extrapolation
%% =========================
%% Correction
%% =========================
%% nonlinear measurement eq
h = hk([uav_init_pos, h0], [uav_actual_pos, h0], [X_k; 0]);
%% Jacobian of nonlinear measurement eq.
H = JH([uav_init_pos, h0],[uav_actual_pos, h0], [X_k; 0]);
%% Equation 3: Innovation
v = alpha_meas - h;
%% ===============================================
%% Equation 4: Innovation
S = H * P_s * H' + R;
%% Equation 5: Kalman gain
K = (P_s * H') / (S);
%% Equation 6: State update
X = X_s + K * (v);
%% Equation 7: Error covariance update
P_f = (eye(size(K * H)) - K * H) * P_s;

%Time update

% X_s = F * X_k;
% P_s = F * P_k * F' + G * Q * G';
% 
% %Measurement update
% Hk = JH([uav_init_pos, h0], [uav_actual_pos, h0], [X_s; 0]);
% 
% K = P_s * Hk' / (Hk * P_s * Hk' + R)
% X = X_s + K * (alpha_meas - hk([uav_init_pos, h0], [uav_actual_pos, h0], [X_s; 0]))
% P_f = (eye(2) - K * Hk) * P_s;

%% ===============================================
%% h(X): Nonlinear measurement eq
function h=hk(uav_init_pos, uav_actual_pos, X_predicted)
h=norm(X_predicted - uav_init_pos')^2 / norm(X_predicted - uav_actual_pos')^2;
%% partial_h/partial_X: Jacobian of measurement eq
function H=JH(uav_init_pos, uav_actual_pos, X_predicted)
x=X_predicted(1);
y=X_predicted(2);
a = norm(X_predicted - uav_init_pos')^2;
b = norm(X_predicted - uav_actual_pos')^2;
H=[(2 * (x - uav_init_pos(1)) * b - 2 * (x - uav_actual_pos(1)) * a) / b^2,...  
   (2 * (y - uav_init_pos(2)) * b - 2 * (y - uav_actual_pos(2)) * a) / b^2];

% function H=JH(s1, sk, xk)
% s1 = s1';
% sk = sk';
% H(1, 1) = 2 / norm(xk - sk) ^ 4 * ((xk(1) - s1(1)) * norm(xk - sk) ^ 2 - (xk(1) - sk(1)) * norm(xk - s1) ^ 2);
% H(1, 2) = 2 / norm(xk - sk) ^ 4 * ((xk(2) - s1(2)) * norm(xk - sk) ^ 2 - (xk(2) - sk(2)) * norm(xk - s1) ^ 2);
