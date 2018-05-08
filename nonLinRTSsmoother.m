function [xs, Ps, xf, Pf, xp, Pp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, sigmaPoints, type)
%NONLINRTSSMOOTHER Filters measurement sequence Y using a
% non-linear Kalman filter.
%
%Input:
%   Y           [m x N] Measurement sequence for times 1,...,N
%   x_0         [n x 1] Prior mean for time 0
%   P_0         [n x n] Prior covariance
%   f                   Motion model function handle
%   T                   Sampling time
%   Q           [n x n] Process noise covariance
%   S           [n x N] Sensor position vector sequence
%   h                   Measurement model function handle
%   R           [n x n] Measurement noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xf          [n x N]     Filtered estimates for times 1,...,N
%   Pf          [n x n x N] Filter error convariance
%   xp          [n x N]     Predicted estimates for times 1,...,N
%   Pp          [n x n x N] Filter error convariance
%   xs          [n x N]     Smoothed estimates for times 1,...,N
%   Ps          [n x n x N] Smoothing error convariance

%% First the entire sequence is filtered

% Parameters
N = size(Y,2);
n = length(x_0);
%m = size(Y,1);

% Data allocation
xp = zeros(n,N);
xf = zeros(n,N);
Pp = zeros(n,n,N);
Pf = zeros(n,n,N);

[xp(:,1), Pp(:,:,1)] = nonLinKFprediction(x_0, P_0, f,T, Q,sigmaPoints, type);
[xf(:,1), Pf(:,:,1)] = nonLinKFupdate(xp(:,1), Pp(:,:,1), Y(:,1),S(:,1), h, R,sigmaPoints, type);

for k=2:N
    [xp(:,k), Pp(:,:,k)] = nonLinKFprediction(xf(:,k-1), Pf(:,:,k-1), f,T, Q,sigmaPoints, type);
    [xf(:,k), Pf(:,:,k)] = nonLinKFupdate(xp(:,k), Pp(:,:,k), Y(:,k),S(:,k), h, R,sigmaPoints, type);
end


%% Smoothing
%no future information is available at time k=N so smoothing is equal to
%filter
xs(:,N) = xf(:,N);
Ps(:,:,N) = Pf(:,:,N);

for k=N-1:-1:1
    
    xs_kplus1= xs(:,k+1);
    Ps_kplus1 = Ps(:,:,k+1);
    xf_k = xf(:,k);
    Pf_k = Pf(:,:,k);
    xp_kplus1 = xp(:,k+1);
    Pp_kplus1 = Pp(:,:,k+1);
    
    [xs(:,k), Ps(:,:,k)] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1,f, T, sigmaPoints, type);
    
    % We have offered you functions that do the non-linear Kalman prediction and update steps.
    % Call the functions using
    % [xPred, PPred] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
    % [xf, Pf] = nonLinKFupdate(xPred, PPred, Y, S, h, R, sigmaPoints, type);
    
end
end