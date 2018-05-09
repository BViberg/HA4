function [Xr, Wr, j] = resampl(X, W)
%RESAMPLE Resample particles and output new particles and weights.
% resampled particles. 
%
%   if old particle vector is x, new particles x_new is computed as x(:,j)
%
% Input:
%   X   [n x N] Particles, each column is a particle.
%   W   [1 x N] Weights, corresponding to the samples
%
% Output:
%   Xr  [n x N] Resampled particles, each corresponding to some particle 
%               from old weights.
%   Wr  [1 x N] New weights for the resampled particles.
%   j   [1 x N] vector of indices refering to vector of old particles

% Generates the segmented numberline from 0 to 1 (upper edge not included since it will never be graer than the edge)
cumulW=cumsum(W);
segment=[0,cumulW(1:end-1)];
N=size(X,2);
samples = rand([N 1]);
%selectionMatrix = samples>cumulW;
j=zeros(1,N);
for i=1:N
    %Check which interval the sample belongs to
    
    % "segment" is strictly monotonicly increasing, hence the index
    % corresponding to the last time the sample is greater than "segment"
    % thats the interval it belongs to
    j(i)=find(segment<samples(i),1,'last');
end
% after resampling all particles have equal weight
Wr=repmat([1/N],[1 N]);
Xr=X(:,j);
end
