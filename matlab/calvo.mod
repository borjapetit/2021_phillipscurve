@#define IRF_periods=36

options_.noprint=1;

var    y pi pi_w w r m mg ;
varexo em;

parameters  theta theta_w gamma elast beta rhom eta epsilon_w ;

@#include "parameters.txt"

gamma     = 2;
elast     = 0.5;
beta      = (1/1.04)^(1/12);
rhom      = 0.8;
eta       = 1/(1/beta-1);       % see Jordi p.27
epsilon_w = 7;


model;

#lambda    = (1-theta)*(1-theta*beta)/theta;                              % Jordi, p. 121
#kappa     = 0*lambda;                                                    % Jordi, p. 126
#lambda_w  = (1-theta_w)*(1-theta_w*beta)/(theta_w*(1+epsilon_w*elast));  % Jordi, p.125
#kappa_w   = lambda_w*(gamma + elast);                                    % Jordi, p.126

y = y(+1) - 1/gamma*(r - pi(+1));
pi = beta*pi(+1) + kappa*y + lambda*w;               % Jordi eq. 15, p. 126
pi_w = beta*pi_w(+1) + kappa_w*y - lambda_w*w;       % Jordi eq. 17, p. 126
w = w(-1) + pi_w - pi;                               % Jordi eq. 18, p. 127
m = gamma*y - eta*r;                                 % Jordi p.27
m - m(-1)= mg - pi;  % real money growth = nominal money gr. less inflation
mg = rhom*mg(-1) + em;

end;

shocks;
var em = 1;
end;

steady;
check;

stoch_simul(order=1,irf=@{IRF_periods},nograph) y pi;

