@#define IRF_periods=36

var    y pi pi_w w r m mg; 
varexo em; 
 
parameters  theta theta_w gamma elast beta rhom eta epsilon_w ;
   
theta     = 0.84846788;  % probability of NOT changing price
theta_w   = 0.705382;  % probability of NOT changing wages
gamma     = 2; 
elast     = 0.5; 
beta      = (1/1.04)^(1/12);
rhom      = 0.8;
eta       = 1/(1/beta-1); % see Jordi p.27
epsilon_w = 7;


model; 
#lambda    = (1-theta)*(1-theta*beta)/theta;                 % Jordi, p. 121
#kappa     = 0*lambda;                                       % Jordi, p. 126    
#lambda_w  = (1-theta_w)*(1-theta_w*beta)/(theta_w*(1+epsilon_w*elast));  % Jordi, p.125
#kappa_w   = lambda_w*(gamma + elast);                                    % Jordi, p.126 

y = y(+1) - 1/gamma*(r - pi(+1));

pi = beta*pi(+1) + kappa*y + lambda*w;               % Jordi eq. 15, p. 126

pi_w = beta*pi_w(+1) + kappa_w*y - lambda_w*w;       % Jordi eq. 17, p. 126

w = w(-1) + pi_w - pi;                               % Jordi eq. 18, p. 127

m = gamma*y - eta*r;                                           % Jordi p.27

m - m(-1)= mg - pi;  % real money growth = nominal money gr. less inflation

mg = rhom*mg(-1) + em;

end; 
 
shocks; 
var em = 1; 
end; 
 
steady;
% check;


stoch_simul(order=1,irf=@{IRF_periods},nograph) y pi;

load V10_irf;

IRF_weighting = eye(2*@{IRF_periods}-1);

IRF_empirical = [irf.consumption(2:@{IRF_periods}+1)' irf.pinflation(3:@{IRF_periods}+2)'];

x_start=[theta theta_w]; 

%make sure Dynare does not print out stuff during runs
options_.nomoments=1;
options_.nofunctions=1;
options_.nograph=1;
options_.verbosity=0;

%set noprint option to suppress error messages within optimizer
options_.noprint=1;

% set csminwel options
H0 = 1e-2*eye(length(x_start)); %Initial Hessian 
crit = 1e-8; %Tolerance
nit = 1000;  %Number of iterations

x_opt_hat = NaN(2,4);

load V10_irf;

IRF_empirical = [irf.consumption(2:@{IRF_periods}+1)' irf.pinflation(3:@{IRF_periods}+2)'];

[fhat,x_opt_hat(:,1)] = csminwel(@IRF_matching_objective,x_start,H0,[],crit,nit,IRF_empirical,IRF_weighting);

load V12_irf;

IRF_empirical = [irf.consumption(2:@{IRF_periods}+1)' irf.pinflation(3:@{IRF_periods}+2)'];

[fhat,x_opt_hat(:,2)] = csminwel(@IRF_matching_objective,x_start,H0,[],crit,nit,IRF_empirical,IRF_weighting);

load V13_irf;

IRF_empirical = [irf.consumption(2:@{IRF_periods}+1)' irf.pinflation(3:@{IRF_periods}+2)'];

[fhat,x_opt_hat(:,3)] = csminwel(@IRF_matching_objective,x_start,H0,[],crit,nit,IRF_empirical,IRF_weighting);

load V14_irf;

IRF_empirical = [irf.consumption(2:@{IRF_periods}+1)' irf.pinflation(3:@{IRF_periods}+2)'];

[fhat,x_opt_hat(:,4)] = csminwel(@IRF_matching_objective,x_start,H0,[],crit,nit,IRF_empirical,IRF_weighting);


%disp(1-x_opt_hat')

