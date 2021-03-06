function [E0,Rs] = DMC(D,dims,V,Rs,E0,dt,M0,steps,bSteps,a)
%DMC Implements the basic diffusion Monte Carlo algorithm.
%   See write-up for detailed explanation and derivation.
M = M0;                           % Current number of walkers (in block)
for block = 1:steps/bSteps        % Repeat for the correct number of blocks
    for step = 1:bSteps           % Repeat for each step in a block
        numWalkers = zeros(M,1);  % Tracks number of walker copies
        for i = 1:M               % For each configuration in Rs
            R = Rs(:,i);          % Retrieve R value
            Rp = R+normrnd(0,sqrt(2*D*dt),dims,1);     % Randomly change R
            weight = exp(-(dt/2)*(V(Rp)+V(R))+dt*E0);  % Branching weight
            numWalkers(i) = floor(weight+rand);        % Number of copies
            Rs(:,i) = Rp;                              % Update R value
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