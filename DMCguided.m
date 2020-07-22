function [E0,Rs] = DMCguided(D,dims,PT,EL,F,E0,dt,M0,steps,bSteps,a)
%DMCGUIDED Implements the guided diffusion Monte Carlo algorithm.
%   See write-up for detailed explanation and derivation.
M = M0;                           % Current number of walkers (in block)
PT2 = @(R) PT(R)^2;               % Trial wave function squared
mult = 10;                        % Multiple for Metropolis sampling
sd = 0.1;                         % St. dev. for Metropolis random step
Rs = zeros(dims,mult*M);          % Allocate array of R values
i = 1;                            % Index to track current R value
while i <= mult*M                 % Repeat until enough Rs are collected
    r = Rs(:,i);                  % Retrieve last accepted R value
    rp = normrnd(r,sd,dims,1);    % Random step
    if PT2(rp)/PT2(r) > rand      % Decide whether to accept or reject
        Rs(:,i+1) = rp;           % If accepted, record R value
        i = i+1;                  % Increment index
    end
end
Rs = Rs(:,(1:M)*mult);            % Select subset of Rs (sampled from PT2)
for block = 1:steps/bSteps        % Repeat for the correct number of blocks
    for step = 1:bSteps           % Repeat for each step in a block
        numWalkers = zeros(M,1);  % Tracks number of walker copies
        for i = 1:M               % For each configuration in Rs
            R = Rs(:,i);          % Retrieve R value
            Rp = R+D*dt*F(R)+normrnd(0,sqrt(2*D*dt),dims,1);  % Step R
            weight = exp(-(dt/2)*(EL(R)+EL(Rp))+dt*E0);  % Branching weight
            numWalkers(i) = floor(weight+rand);          % Number of copies
            Rs(:,i) = Rp;                                % Update R value
        end
        newM = sum(numWalkers);           % New number of walkers
        newRs = zeros(dims,newM);         % Allocate new list of walkers
        newi = 1;                         % Current index in newRs
        for i = 1:M                       % For each old walker R
            for j = 1:numWalkers(i)       % Repeat for each copy of R
                newRs(:,newi) = Rs(:,i);  % Save the copy as a newR
                newi = newi+1;            % Increment index in newRs
            end
        end
        Rs = newRs;                       % Update current Rs list
        M = newM;                         % Update number of walkers
    end
    E0 = E0-a*(M-M0)/(dt*M0);             % Adjust estimate for E0
    M0 = M;                               % Update intial number of walkers
end
end