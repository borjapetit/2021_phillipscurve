
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE TAKES THE JACOBIAN COMPUTED WITH FORTRAN, AND COMPUTE THE IRFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

% Working directory
cd ' %% Your Own Directory '
mypath = [pwd '/textfiles/'] ;

% Define the versions you want to solve
cases = [17 18];

% Write jacobian in .mat file
fprintf('\n CONVERTING TEXT FILES INTO MAT FILES \n\n');
for version_kappas = cases
    fprintf(' Converting version %d\n',version_kappas);
    filetxt = [mypath '_dyn/V' num2str(version_kappas) '_dyn.txt'] ;
    filemat = [mypath '_dyn/V' num2str(version_kappas) '_dyn.mat'] ;
    if isfile(filetxt)
        sol_dyn = importdata(filetxt,' ',0);
        if isfile(filemat) ; delete(filemat) ; end
        eval(['save ' filemat ' sol_dyn']);
        delete(filetxt)
    end
end

% Compute the IRFs
for version_kappas = cases
    solde_dyn_version(version_kappas,mypath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RES = solde_dyn_version(version_kappas,mypath)

close all ; tic ;

fprintf('\n DYNAMICS OF VERSION %d - STARTED AT %s\n\n',version_kappas,datetime('now','Format','HH:mm'));

% Figure size parameters
fig_w = 10.0 ;
fig_h = 4.00 ;

% Program execution parameters
tol      = sqrt(eps) ;
eqcutoff = tol*1000  ;
jacstep  = tol       ;

% Model parameters
run('parameters.m');

% Define inflation rate
if version_kappas>6
    mu0 = mu ; version_inf = version_kappas - 10*floor(version_kappas/10);
    if version_inf==1 ; mu = (mu0^12 - 0.03)^(1/12) ; end
    if version_inf==2 ; mu = (mu0^12 - 0.02)^(1/12) ; end
    if version_inf==3 ; mu = (mu0^12 + 0.02)^(1/12) ; end
    if version_inf==4 ; mu = (mu0^12 + 0.06)^(1/12) ; end
    if version_inf==5 ; mu = (mu0^12 - 0.04)^(1/12) ; end
    if version_inf==6 ; mu = (mu0^12 - 0.01)^(1/12) ; end
    if version_inf==7 ; mu =        1.046342^(1/12) ; end
    if version_inf==8 ; mu =        1.020054^(1/12) ; end
end

% Loading solution from Fortran
fprintf(' Loading solution from Fortran...\n');

% Extract steady-state solution
sol_ss = importdata([mypath '_ss/V' num2str(version_kappas) '_ss.txt'],' ',0);
run('extract_ss.m');
clear sol_ss

% Extract histograms of price and wage changes (from data)
pricehistdatap = importdata([pwd 'textfiles/data_pdfprices.txt'],' ',0) ;
pricehistdataw = importdata([pwd 'textfiles/data_pdfwages.txt'],' ',0) ;

% Extract dynamic solution
eval(['load ' mypath '_dyn/V' num2str(version_kappas) '_dyn sol_dyn']);
run('extract_dyn.m');
clear sol_dyn

% Plot histogram of price/wage changes
fprintf(' Printing histograms... \n');

fig_hist = figure;
    ip = 1:1:(2*nump-1) ; ip = ip' ; xp = [ip';ip']; yp = [pricehistdatap';pricehistdatap'];
    iw = 1:1:(2*numw-1) ; iw = iw' ; xw = [iw';iw']; yw = [pricehistdataw';pricehistdataw'];
subplot(1,2,1)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,hist_p,'b','LineWidth',2);
    title("Prices")
    grid on
    legend('  Data', '  Model')
    hold off
subplot(1,2,2)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,hist_w,'b','LineWidth',2);
    title("Wages")
    grid on
    hold off
set(fig_hist,'PaperSize',[2*fig_w 2*fig_h],'PaperPosition',[0 0 2*fig_w 2*fig_h])

% Save figure
filename = [mypath '_ss/_figs_hist/V' num2str(version_kappas) '_hist.pdf'];
if isfile(filename); delete(filename); end
print(fig_hist,filename,'-dpdf')
clear fig_hist filename ip iw xp xw yp yw
close all

% Compute Kelin's QZ decomposition
fprintf(' Start Klein solution... \n');
[JUMPS,STATEDYNAMICS,stableeigs,unstableeigs] = ...
    kleinsolve(ajac,bjac,cjac,djac,phiMAT,nPdist+nWdist+nss,eqcutoff);
tocdyn = toc;
fprintf(' Elapsed time solving dynamics: %d min and %d secs \n\n',...
    floor(tocdyn/60),floor(tocdyn-60*floor(tocdyn/60)));
clear tocjacobian tockelin tocdyn


% Set parameters for IRF computation
clear STATEHISTORY; clear JumpHistory; TT = 200;

% Simulate monetary shock
Shocks      = zeros(nz+nPdist+nWdist+nss,TT);
Shocks(1,2) = jacstep;

% Solve for the variable paths
States      = zeros(nz+nPdist+nWdist+nss,TT);
Jumps       = zeros(nV+nL+nsj,TT);
for time = 2:TT
    States(:,time) = STATEDYNAMICS*States(:,time-1) + Shocks(:,time);
    Jumps(:,time)  = JUMPS*States(:,time);
end

% Shoks paths
R_path     = States(1,:);
S_path     = States(2,:);
Z_path     = States(3,:);

% State variables paths
Pdist_path = Pdist(:)*ones(1,TT) + States(nz+1:nz+nPdist,:);
Wdist_path = Wdist(:)*ones(1,TT) + States(nz+nPdist+1:nz+nPdist+nWdist,:);
mlag_path  = mbar                + States(nz+nPdist+nWdist+1,:);
w_path     = wbar                + States(nz+nPdist+nWdist+2,:);
Inf_path   = mu                  + States(nz+nPdist+nWdist+3,:);

% Jump variables paths
V_path     = V(:)*ones(1,TT)     + Jumps(1:nV,:);
L_path     = L(:)*ones(1,TT)     + Jumps(nV+1:nV+nL,:);
C_path     = Cbar                + Jumps(nV+nL+1,:);
N_path     = Nbar                + Jumps(nV+nL+2,:);
Inf_w_path = [mu,w_path(1,2:TT)./w_path(1,1:TT-1).*Inf_path(1,2:TT)];

% Compute IRFs
shocksize       = R_path(2)                             ;
irf.money       = R_path/shocksize                      ;
irf.pinflation  = (Inf_path   - mu   )/shocksize        ;
irf.consumption = (C_path     - Cbar )/(Cbar*shocksize) ;
irf.labor       = (N_path     - Nbar )/(Nbar*shocksize) ;
irf.winflation  = (Inf_w_path - mu   )/shocksize        ;
irf.wage        = (w_path     - wbar )/(wbar*shocksize) ;

% Save dynamic results
if isfile([mypath '_irfs/V' num2str(version_kappas) '_irf.mat'])
    eval(['delete ' mypath '_irfs/V' num2str(version_kappas) '_irf.mat']);
end
eval( ['save ' mypath '_irfs/V' num2str(version_kappas) '_irf irf STATEDYNAMICS JUMPS '] );

% Plot IRFs
fig_irf = figure;
hold on ; linescolor = '-o' ; plot_irf
set(fig_irf,'PaperSize',[3*fig_w 3*fig_h],'PaperPosition',[0 0 3*fig_w 3*fig_h])

% Save IFR figure
filename = [mypath '_irfs/_figs_irfs/V' num2str(version_kappas) '_irf.pdf'];
if isfile(filename) ; delete(filename); end
print(fig_irf,filename,'-dpdf')
clear fig_irf filename
close all

% Finished
RES = irf.wage(end) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
