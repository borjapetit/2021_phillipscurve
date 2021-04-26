
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix and vector sizes

nump   = 36 ; numw = 71 ;   % # of points in the grid of prices/wages
nums   = 25 ; numz = 51 ;   % # of points in the grid of firm's/worker's productivity

nsj    = 2  ;               % # of scalar jump variables
nss    = 3  ;               % # of scalar state variables
nz     = 3  ;               % # of exogenous shock processes

nPdist = nump*nums;         % # of points of the distribution
nWdist = numw*numz;         % # of points of the distribution
nV     = nump*nums;         % # of points of the value function
nL     = numw*numz;         % # of points of the value function

Ntot   = 2*nPdist+2*nWdist+nss+nsj; % # of variables (excluding shocks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters

delta       = 1-1.04^(-1/12);        % Monthly time discount factor
gamma       = 2;                     % CRRA exponent
nu          = 1;                     % Money demand parameter
mu          = (1+4.2604e-004)^4;     % Monthly rate of money growth (~2%/yr)
Rss         = mu/(1-delta);          % SS interest rate

rho_r       = 0.800000000000000;     % Monthly persistence of money shock 
rho_z       = 0.970000028289202;     % Monthly persistence of firms TFP
rho_s       = 0.644109165109414;     % Monthly persistence of workers TFP

phiMAT      = zeros(nz,nz) ;         % Matrix of shock persistence 
phiMAT(1,1) = rho_r ;                % Persistence of money growth shock
phiMAT(2,2) = rho_s ;                % Persistence of productivty shocks
phiMAT(3,3) = rho_z ;                % Persistence of workers productivity shocks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price and age grids:

p_grid = [ -0.350 : 0.70/(nump-1) : 0.35 ]' ;
w_grid = [ -0.350 : 0.70/(numw-1) : 0.35 ]' ;



