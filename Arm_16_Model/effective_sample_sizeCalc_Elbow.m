clear
clc

% calculate convergence diagnostics Rhat and effective sample size
% https://mc-stan.org/docs/2_26/reference-manual/effective-sample-size-section.html

% Author: Russell T. Johnson, rtjohnso@usc.edu
% Last Edited: 6-18-21

% load chain from file
load  chain_results_20220115T023505

[n_iter, param, M] = size(chain); 
burn_in = n_iter*0.5; % specify the burn in 


for i = 1:M 
   chain_reduce(:,:,i) = chain(burn_in+1:100:end,:,i);
end


% N = n_iter-burn_in;
N = size(chain_reduce,1); 

% calculate R Hat
for p = 1:param

    % chain_corr = chain(:,:,choose_result); 
    theta = chain_reduce(:,p,:); 

    % between chain variance
    thetadotm = (1/N)*sum(theta,1);
    thetadotdot = (1/M)*sum(thetadotm); 
    samp(p) = thetadotdot;
    for i = 1:M
        diff(i) = (thetadotm(i) - thetadotdot)^2;
    end
    B = (1/(M-1))* (sum(diff));

    % within chain variance
    for i = 1:M
        for j = 1:N
            temp(j,i) = (theta(j,:,i)-thetadotm(i))^2;
        end
        smsquare(i) = (1/(N-1))*(sum(temp(j,i)));
    end
    W = (1/M) * sum(smsquare);

    % variance estimator
    varhat = ((N-1)/N)*W + (1/N)*B; 

    % potential scale reduction statistic 
    Rhat(p) = sqrt(varhat/W); 

end


%% from Stan documentation
for p = 1:param
    % choose_param = 1; 
    theta = chain_reduce(:,p,:); 

    thetadotm = (1/N)*sum(theta,1);
    thetadotdot = (1/M)*sum(thetadotm); 
    for i = 1:M
        diff(i) = (thetadotm(i) - thetadotdot)^2;
    end
    B = (1/(M-1))* (sum(diff));

    % within chain variance
    for i = 1:M
        for j = 1:N
            temp(j,i) = (theta(j,:,i)-thetadotm(i))^2;
        end
        smsquare(i) = (1/(N-1))*(sum(temp(j,i)));
    end
    W = (1/M) * sum(smsquare);

    % variance estimator
    varhat = ((N-1)/N)*W + (1/N)*B; 

    for i = 1:M
         AutoCorr_msd(:,i) = autocorr(theta(:,:,i),N-1); 
    end
    % PTM = sum(AutoCorr_msd,2);

    for i = 1:M
        quot(:,i) = smsquare(i).*AutoCorr_msd(:,i);
    end

    for i = 1:N
        PThat(i) = 1 - ((W - (1/M)*sum(quot(i,:)))/varhat);
    end

    acf_sum = 10; 
    t = 1; 
    while (acf_sum >= 0)
       t = t+1;
       acf_sum = AutoCorr_msd(t) + AutoCorr_msd(t+1); 
    end

    for k = 1:t
        autocorr_sum(k) = AutoCorr_msd(k) + AutoCorr_msd(k+1); 
    end

    That = -1 + 2 * sum(autocorr_sum); 

    NEFF(p) = (M*N)/That; 

end








