
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 4 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

% Working directory
cd ' %% Your Own Directory '

% Size parameters
run('matlab/parameters.m');

% Load data on price and wage changes
pricehistdatap = importdata([pwd '/textfiles/data_pdfprices.txt'],' ',0) ;
pricehistdataw = importdata([pwd '/textfiles/data_pdfwages.txt'],' ',0) ;

% Create plot for data
ip = -0.70:1.40/(2*nump-2):0.70 ; ip = ip' ; xp = [ip';ip']; yp = [pricehistdatap';pricehistdatap'];
iw = -0.70:1.40/(2*numw-2):0.70 ; iw = iw' ; xw = [iw';iw']; yw = [pricehistdataw';pricehistdataw'];

% Load model solutions for baseline, FP and FW versions
sol_ss1 = importdata([pwd '/textfiles/_ss/V10_ss.txt'],' ',0);
sol_ss3 = importdata([pwd '/textfiles/_ss/V30_ss.txt'],' ',0);
sol_ss5 = importdata([pwd '/textfiles/_ss/V50_ss.txt'],' ',0);

% Extract model-generated histograms for versions 2, 3, 4 and 5
ph1 = sol_ss1(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);
ph3 = sol_ss3(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);
ph5 = sol_ss5(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);

% Create the figure
fig_hist = figure;
subplot(1,2,1)
    hold on
    %area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,ph1(1:nump*2-1),'r','LineWidth',1);
    stairs(ip,ph3(1:nump*2-1),'b','LineWidth',1);
    title("Flexible prices",'FontWeight','normal')
    xlabel("Log price changes")
    legend(' Base.',' FP')
    grid on
    box on
    xlim([-0.5 0.5])
    xticks([-0.5 -0.25 0 0.25 0.5])
    hold off
subplot(1,2,2)
    hold on
    %area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,ph1(nump*2:nump*2+numw*2-2),'r','LineWidth',1);
    stairs(iw,ph5(nump*2:nump*2+numw*2-2),'b','LineWidth',1);
    title("Flexible wages",'FontWeight','normal')
    xlabel("Log wage changes")
    legend(' Base.',' FW')
    grid on
    box on
    xlim([-0.25 0.25])
    hold off

set(fig_hist,'PaperSize',[20 7],'PaperPosition',[0 0 20 7])
savefig('figures/figs/fig_4.fig')
print(fig_hist,[pwd '/figures/pdfs/fig_4.pdf'],'-dpdf')
