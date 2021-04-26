
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 3 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

% Working directory
cd '/Users/borjapetit/Dropbox/projects/2016_lpw/code_public/codes/'

sol_ss = importdata([pwd '/textfiles/_ss/V10_ss.txt'],' ',0);
lambda = importdata([pwd '/textfiles/_ss/_V10_lambda_ss.txt'],' ',0);
Pi     = importdata([pwd '/textfiles/_ss/_V10_pi_ss.txt'],' ',0);
rho    = importdata([pwd '/textfiles/_ss/_V10_rho_ss.txt'],' ',0);
Pi_w   = importdata([pwd '/textfiles/_ss/_V10_piw_ss.txt'],' ',0);

run('matlab/parameters.m');
run('matlab/extract_ss.m');

Pi      = reshape(Pi       ,nums,nump) ;
lambda  = reshape(lambda   ,nums,nump) ;
Pi_w    = reshape(Pi_w,numw,numz,numw) ;
rho     = reshape(rho      ,numz,numw) ;

auxv   = cumsum(sum(Pdist)');
[~,a1] = min(abs(auxv-0.25*ones(25,1)));
[~,a2] = min(abs(auxv-0.50*ones(25,1)));
[~,a3] = min(abs(auxv-0.75*ones(25,1)));

auxv    = cumsum(sum(Wdist)');
[i1,b1] = min(abs(auxv-0.25*ones(51,1)));
[i2,b2] = min(abs(auxv-0.50*ones(51,1)));
[i3,b3] = min(abs(auxv-0.75*ones(51,1)));

% Create the figure
fig = figure;
subplot(2,2,1)
    hold on
    plot(p_grid,Pi(a1,:),'LineWidth',1.5)
    plot(p_grid,Pi(a2,:),'LineWidth',1.5)
    plot(p_grid,Pi(a3,:),'LineWidth',1.5)
    title("\rm Panel A. Firms: Logit probabilities")
    xlabel("(log) Prices")
    xlim([-0.2 0.2])
    legend(' p25',' p50',' p75','Location','NorthWest')
    grid on
    box on
    hold off
subplot(2,2,2)
    hold on
    plot(p_grid,lambda(a1,:),'LineWidth',1.5)
    plot(p_grid,lambda(a2,:),'LineWidth',1.5)
    plot(p_grid,lambda(a3,:),'LineWidth',1.5)
    title("\rm Panel B. Firms: Probabilities of adjustment")
    xlabel("(log) Prices")
    ylim([0 1])
    xlim([-0.2 0.2])
    grid on
    box on
    hold off
subplot(2,2,3)
    hold on
    plot(w_grid,squeeze(Pi_w(:,b1,35)),'LineWidth',1.5)
    plot(w_grid,squeeze(Pi_w(:,b2,35)),'LineWidth',1.5)
    plot(w_grid,squeeze(Pi_w(:,b3,35)),'LineWidth',1.5)
    title("\rm Panel C. Workers: Logit probabilities")
    xlabel("(log) Wages")
    xlim([-0.1 0.3])
    grid on
    box on
    hold off
subplot(2,2,4)
    hold on
    plot(w_grid,rho(b1,:),'LineWidth',1.5)
    plot(w_grid,rho(b2,:),'LineWidth',1.5)
    plot(w_grid,rho(b3,:),'LineWidth',1.5)
    title("\rm Panel D. Workers: Probabilities of adjustment")
    xlabel("(log) Wages")
    ylim([0 1])
    xlim([-0.1 0.3])
    grid on
    box on
    hold off

set(fig,'PaperSize',[20 15],'PaperPosition',[0 0 20 15])
savefig([pwd '/figures/figs/fig_3.fig'])
print(fig,[pwd '/figures/pdfs/fig_3.pdf'],'-dpdf')
