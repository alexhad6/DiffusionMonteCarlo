function [E0,R] = DMC(D,dims,pT,EL,F,E0,dt,M0,steps,bSteps,a)
% Sample initial configuration from pT using Metropolis algorithm
mult = 10;
sd = 0.1;
R = sampleFrom(@(r)pT(r)^2,dims,M0,mult,sd);

% Run DMC simulation
M = M0;                           % Current number of walkers (in block)
for block = 1:steps/bSteps
    for step = 1:bSteps
        numWalkers = zeros(M,1);  % Tracks number of walker
        
        % Step walkers forward
        for i = 1:M
            r = R(:,i);
            rp = r+D*dt*F(r)+normrnd(0,sqrt(2*D*dt),dims,1);
            weight = exp(-(dt/2)*(EL(r)+EL(rp))+dt*E0);
            numWalkers(i) = floor(weight+rand);
        end
        
        % Update number of walkers
        newM = sum(numWalkers);
        newR = zeros(dims,sum(newM));
        newi = 1;
        for i = 1:M
            for j = 1:numWalkers(i)
                newR(:,newi) = R(:,i);
                newi = newi+1;
            end
        end
        R = newR;
        M = size(R,2);
    end
    E0 = E0-a*(M-M0)/(dt*M0);
    M0 = M;
end
end