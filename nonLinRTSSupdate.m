function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
    Ps_kplus1, ...
    xf_k, ...
    Pf_k, ...
    xp_kplus1, ...
    Pp_kplus1, ...
    f, ...
    T, ...
    sigmaPoints, ...
    type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

switch type
    case 'EKF'
        [fxhat,Fxhat]=f(xf_k);
        G = Pf_k*Fxhat'/Pp_kplus1; 
        xs = xf_k+G*(xp_kplus1-fxhat);
        Ps = Pf_k-G*(Pp_kplus1-Ps_kplus1)*G';
        
    case 'UKF'
        
    case 'CKF'
            
    
    otherwise
        error('Invalid type')
end
        
end