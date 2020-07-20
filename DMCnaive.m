function [E0,R] = DMCnaive(D,dims,V,R,E0,dt,M0,steps,bSteps,a)
M = M0;                       % Current number of walkers (in block)
for block = 1:steps/bSteps
    for step = 1:bSteps
        numWalkers = zeros(M,1);  % Tracks number of walker
        
        % Step walkers forward
        for i = 1:M
            r = R(i,:);
            rp = r+normrnd(0,sqrt(2*D*dt),1,dims);
            weight = exp(-(dt/2)*(V(rp)+V(r))+dt*E0);
            numWalkers(i) = floor(weight+rand);
            R(i,:) = rp;
        end
        
        % Update number of walkers
        newM = sum(numWalkers);
        if newM ~= M
            newR = zeros(sum(newM),dims);
            newi = 1;
            for i = 1:M
                for j = 1:numWalkers(i)
                    newR(newi,:) = R(i,:);
                    newi = newi+1;
                end
            end
            R = newR;
            M = size(R,1);
        end
    end
    E0 = E0-a*(M-M0)/(dt*M0);
    M0 = M;
end
end