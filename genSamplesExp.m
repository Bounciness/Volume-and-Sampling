function [ points ] = genSamplesExp( P,lambda, numSteps, numSamples, warmup )
%GENSAMPLES Compute samples from a polytope, taking numSteps steps of the
%Markov chain before saving a sample point

if nargin < 5
   warmup = 1e4; 
end

if nargin < 4
   numSamples = 1e3; 
end

if nargin < 3
    numSteps = 1e2;
end

K = ConvexBody(P,[],.2,'');
x = zeros(K.dim,1);
resetSlacks(K,x);
points = zeros(size(P.N,1),numSamples);

for i=1:warmup
   x = getNextPoint_Exp(K,x,lambda,1); 
end

h = waitbar(0,'Computing samples...');
for i=1:numSamples
    for j=1:numSteps
        x = getNextPoint_Exp(K,x,lambda,1);
    end
    points(:,i) = P.N*x+P.p_shift;
    waitbar(i/numSamples);
end
close(h);
end