function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

    switch type        
        case 'UKF'
            n=length(x);
            W=zeros(1,2*n+1);
            W(1)=1-n/3; 
            W(2:end)=(1-W(1))/n;
            SP=zeros(n,2*n+1);
            SP(:,1)=x;
            sqrtP=sqrtm(0.5*(P+P'));
            for i=1:n
                SP(:,i+1)=SP(:,1)+sqrt(n/(1-W(1)))*sqrtP(:,i);                
                SP(:,n+i+1)=SP(:,1)-sqrt(n/(1-W(1)))*sqrtP(:,i);
            end
        case 'CKF'
            n=length(x);
            W=zeros(1,2*n)+1/2/n;
            
            SP=zeros(n,2*n);
            sqrtP=sqrtm(0.5*(P+P'));
            for i=1:n
                SP(:,i)=x+sqrt(n)*sqrtP(:,i);                
                SP(:,n+i)=x-sqrt(n)*sqrtP(:,i);
            end
        otherwise
            error('Incorrect type of sigma point')
    end

end