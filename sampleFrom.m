function R = sampleFrom(f,dims,M,mult,sd)
%SAMPLEFROM Samples from the function f using Metropolis algorithm
%   Samples from the function f, which is a wave function squared, so it
%   takes in a dims dimensional vector and outputs a scalar.
multR = zeros(dims,mult*M);
i = 1;
while i <= mult*M
    r = multR(:,i);
    rp = normrnd(r,sd,dims,1);
    if f(rp)/f(r) > rand
        multR(:,i+1) = rp;
        i = i+1;
    end
end
R = multR(:,(1:M)*mult);
end

