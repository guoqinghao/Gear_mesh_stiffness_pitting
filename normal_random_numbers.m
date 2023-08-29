mu = [3.3333+15 0.2682 ];
sigma = [2.4445 0; 0 0.1903];
rng('default')  % For reproducibility
R = mvnrnd(mu,sigma,30);
RR=[R(:,1),-R(:,2)];
plot(RR(:,1),RR(:,2),'*');