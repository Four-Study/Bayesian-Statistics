%% Problem 5 part 2
trueMu = 200;
trueSigma = sqrt(2);
x = trueSigma * randn(100, 1) + trueMu;
[mu,phi] = MH_Normal(x, trueMu, trueSigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  M-H sampling example
%    Implements a Metropolis-Hastings
%     sampler for the 
%     mean and variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [mu,phi] = MH_Normal(x, trueMu, trueSigma)
 
%set basical functionality controls
trace=1;
 
thin = 1;
step = 1; %proposal variance
iterations = 1000000;
burnin = 4000;
 
muInit=0;%trueMu;
sr_phiInit=2;%1/trueSigma;
phiInit=sr_phiInit^2;%1/trueSigma^2;
 
%Note: priors...
%priorMu = 1;
%priorVar = 1/sigma^2; 
 
%%Create space for the stored variables
mu = zeros(1,iterations);
phi = zeros(1,iterations);
sr_phi = zeros(1,iterations);
 
%initialize
mu(1) = muInit;
sr_phi(1) = sr_phiInit;
phi(1) = phiInit;
N = length(x);
 
%basic plotting of the "truth"... true inputs are only used here
figure(1)
clf
plot(trueMu, 1/trueSigma^2,'r.')
xlabel('\mu','fontsize',14)
ylabel('\phi','fontsize',14)
pause
hold on

tic
for i=2:iterations
 
    %%update the mean (Sample from Normal(\bar x),1/(N\phi)), but in an MH
    %%step
    proposalMu = mu(i-1) + randn*step*5;
    proposalPhi_sr = sr_phi(i-1) + randn*step/10;
    proposalPhi = proposalPhi_sr^2;
    %log scale??    
    %proposed = log(prod((1/sqrt(2*pi)*sqrt(phi(i-1))*exp(-phi(i-1)/2*(x-proposalMu).^2))));
    %currentlyAt = log(prod((1/sqrt(2*pi)*sqrt(phi(i-1))*exp(-phi(i-1)/2*(x-mu(i-1)).^2))));
    proposed    = (N/2 - 1)*log(proposalPhi) - proposalPhi * sum((x - proposalMu).^2)/2;
    currentlyAt = (N/2 - 1)*log(phi(i-1)) - phi(i-1) * sum((x - mu(i-1)).^2)/2;
    
    ratio = proposed-currentlyAt;

    if(log(rand)<ratio)
        mu(i) = proposalMu;
        sr_phi(i) = proposalPhi_sr;
        phi(i) = proposalPhi;
    else
        mu(i) = mu(i-1);
        sr_phi(i) = sr_phi(i-1);
        phi(i) = phi(i-1);
    end
 
end
toc
end    

%% Problem 5 part 3
x2 = trueSigma * randn(100, 1) + trueMu;
% lecture code (not attached)
[mu2,phi2] = Gibbs_Normal(x2, trueMu, trueSigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Gibbs sampling example
%    Implements a Gibbs
%     sampler for the 
%     mean and variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [mu,phi] = Gibbs_Normal(x, trueMu, trueSigma)
 
%set basical functionality controls
trace=1;
 
thin = 1;
step = 1; %proposal variance
iterations = 100000;
burnin = 5;
 
muInit= 0;%trueMu;
phiInit=5;%1/trueSigma^2;
 
%Note: priors...
%priorMu = 1;
%priorVar = 1/sigma^2; 
 
 
%%Create space for the stored variables
mu = zeros(1,iterations);
phi = zeros(1,iterations);
 
%initialize
mu(1) = muInit;
phi(1) = phiInit;
N = length(x);
 
%basic plotting of the "truth"... true inputs are only used here
figure(1)
clf
plot(trueMu, 1/trueSigma^2,'r.')
xlabel('\mu','fontsize',14)
ylabel('\phi','fontsize',14)
pause
hold on
 
 
%Gibbs Sample (This really is using an embedded M-H sampler for
%instructional purposes... More on that LATER:)).
tic
for i=2:iterations
 
    %%update the mean (Sample from Normal(\bar x),1/(N\phi)), but in an MH
    %%step
    mu(i) = randn * 1/sqrt(N * phi(i-1)) + mean(x);
    phi(i) = gamrnd( N / 2, 2 / sum((x - mu(i)).^2) );
 
    plotStep=10000; 
    if(mod(i,plotStep)==0)
        figure(1)
        clf
        xlabel('\mu','fontsize',18)
        ylabel('\phi','fontsize',18)
        i
        hold on
        plot(mu(1:i),phi(1:i),'-','color',[100/255,100/255,100/255])
        plot(mu(i-plotStep+1:i),phi(i-plotStep+1:i),'k-')
        plot(trueMu, 1/trueSigma^2,'r*')
        pause(.1)
    end
 
end
toc
    
end    
