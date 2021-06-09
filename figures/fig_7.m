
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 7 and TABLE 8 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic ; format long ; savepwd = pwd;

% Working directory
cd '..'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTUAL DATA

dat      = load([pwd '/textfiles/data_pc.mat']);
ypre_d   = dat.ypre_d   ;
ypost_d  = dat.ypost_d  ;
pipre_d  = dat.pipre_d  ;
pipost_d = dat.pipost_d ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATED DATA

[pipre,ypre,ppre]    = getseries("pre");
[pipost,ypost,ppost] = getseries("post");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE DATA

% zero mean
ypre_d   = ypre_d   - mean(ypre_d);
ypost_d  = ypost_d  - mean(ypost_d);
pipre_d  = pipre_d  - mean(pipre_d);
pipost_d = pipost_d - mean(pipost_d);
ypre     = ypre     - mean(ypre);
ypost    = ypost    - mean(ypost);
pipre    = pipre    - mean(pipre);
pipost   = pipost   - mean(pipost);

% scaling factor for model dispersion
stdpre  =  std(pipre_d)/std(pipre);
stdpost = std(pipost_d)/std(pipost);

% scale model-generated dispersion
ypre   =   ypre.*stdpre  ;
pipre  =  pipre.*stdpre  ;
ypost  =  ypost.*stdpost ;
pipost = pipost.*stdpost ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEST FIT TO SCATTER PLOTS

yfit = [-0.8:0.1:0.8]';

% Estimation - Data
Fitpre_d  = polyfit(ypre_d,pipre_d,1)   ; pifitpre_d  = Fitpre_d(2)  + Fitpre_d(1)*yfit;
Fitpost_d = polyfit(ypost_d,pipost_d,1) ; pifitpost_d = Fitpost_d(2) + Fitpost_d(1)*yfit;

% Estimation - Model
Fitpre  = polyfit(ypre,pipre,1)   ; pifitpre  = Fitpre(2)  + Fitpre(1)*yfit;
Fitpost = polyfit(ypost,pipost,1) ; pifitpost = Fitpost(2) + Fitpost(1)*yfit;

diary 'tables/table_8.txt'

fprintf('\n')
fprintf('\n')
fprintf('    TABLE 8. SLOPE OF THE PHILLIPS CURBE. DATA VS. MODEL')
fprintf('\n')
fprintf('\n')
fprintf('    ***************************************************** \n')
fprintf('    PC Slope   1980-2000   2000-2020     Diff    %% Change \n')
fprintf('    ***************************************************** \n')
fprintf('    Data        %6.4f      %6.4f     %6.4f    %4.2f%% \n', Fitpre_d(1) , Fitpost_d(1) , Fitpost_d(1)-Fitpre_d(1) , 100*Fitpost_d(1)/Fitpre_d(1)-100 )
fprintf('    Model       %6.4f      %6.4f     %6.4f    %4.2f%% \n', Fitpre(1)   , Fitpost(1)   , Fitpost(1)-Fitpre(1)     , 100*Fitpost(1)/Fitpre(1)-100     )
fprintf('    ***************************************************** \n')
fprintf('                    %% observed decline explained = %4.2f%% \n', 100*(Fitpost(1)/Fitpre(1)-1)/(Fitpost_d(1)/Fitpre_d(1)-1) )
fprintf('      \n')
fprintf('      \n')

diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz     = 15;
tickss = [-0.08 -0.06 -0.04 -0.02 0.00 0.02 0.04 0.06 0.08];

f1 = figure; hold on;
subplot(2, 2, 1);hold on;
    title('Phillips Curve 1980-2000. Data','FontWeight','normal')
    scatter(ypre_d,pipre_d,sz,'b')
    plot(yfit,pifitpre_d,'b','LineWidth',1.0)
    xticks(tickss); xlim([-0.06 0.06]) ; xlabel('Output gap')
    yticks(tickss); ylim([-0.06 0.06]) ; ylabel('Change in inflation')
    grid on ; box on
subplot(2, 2, 2);hold on;
    title('Phillips Curve 2000-2020. Data','FontWeight','normal')
    scatter(ypost_d,pipost_d,sz,'b')
    plot(yfit,pifitpost_d,'b','LineWidth',1.0)
    xticks(tickss); xlim([-0.06 0.06]) ; xlabel('Output gap')
    yticks(tickss); ylim([-0.06 0.06]) ; ylabel('Change in inflation')
    grid on ; box on
subplot(2, 2, 3);hold on;
    title('Phillips Curve 1980-2000. Model','FontWeight','normal')
    scatter(ypre,pipre,sz,'r')
    plot(yfit,pifitpre,'r','LineWidth',1.0)
    xticks(tickss); xlim([-0.06 0.06]) ; xlabel('Output gap')
    yticks(tickss); ylim([-0.06 0.06]) ; ylabel('Change in inflation')
    grid on ; box on
subplot(2, 2, 4);hold on;
    title('Phillips Curve 2000-2020. Model','FontWeight','normal')
    scatter(ypost,pipost,sz,'r')
    plot(yfit,pifitpost,'r','LineWidth',1.0)
    xticks(tickss); xlim([-0.06 0.06]) ; xlabel('Output gap')
    yticks(tickss); ylim([-0.06 0.06]) ; ylabel('Change in inflation')
    grid on ; box on
hold off
set(f1,'PaperSize',[20 14],'PaperPosition',[0 0 20 14])
savefig('figures/figs/fig_7.fig')
print(f1,'figures/pdfs/fig_7','-dpdf')

toc ; cd(savepwd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [inf,outp,prices] = getseries(mod)

tol      = sqrt(eps) ;
eqcutoff = tol*1000  ;
jacstep  = tol       ;

run([pwd '/matlab/parameters.m']);

if mod=="pre"
    sol_ss = importdata([pwd '/textfiles/_ss/V17_ss.txt'],' ',0) ;
    load([pwd '/textfiles/_irfs/V17_irf.mat']) ;
    mu = 1.046342^(1/12) ;
end
if mod=="post"
    sol_ss = importdata([pwd '/textfiles/_ss/V18_ss.txt'],' ',0) ;
    load([pwd '/textfiles/_irfs/V18_irf.mat']) ;
    mu = 1.020054^(1/12) ;
end
run([pwd '/matlab/extract_ss.m']);

rng(1);

TT = 1000;

Shocks      = zeros(nz+nPdist+nWdist+nss,TT);
Shocks(1,:) = jacstep*randn(1,TT);
States      = zeros(nz+nPdist+nWdist+nss,TT);
Jumps       = zeros(nV+nL+nsj,TT);
for time = 2:TT
    States(:,time) = STATEDYNAMICS*States(:,time-1) + Shocks(:,time);
    Jumps(:,time)  = JUMPS*States(:,time);
end
Inf_path  = mu   + States(nz+nPdist+nWdist+3,:);
C_path    = Cbar + Jumps(nV+nL+1,:);

inflation = (Inf_path' - mu   ) ;
outputgap = (C_path'   - Cbar ) / ( Cbar ) ;
prices    = cumsum(inflation) ;

j    = 0;
vec1 = NaN(TT,1);
vec2 = NaN(TT,1);
for i=4:3:size(inflation,1) ; j = j + 1 ;
   vec1(j) = exp(prices(i,1))/exp(prices(i-3,1)) - 1.00 ;
   vec2(j) = outputgap(i,1);
end
inf  = vec1(2:j)-vec1(1:j-1) ;
outp = vec2(2:j) ;

if size(inf,1)>80
    inf = inf(end-79:end);
end
if size(outp,1)>80
    outp = outp(end-79:end);
end

prices = prices(end-79:end);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
