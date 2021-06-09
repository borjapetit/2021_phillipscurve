
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 2 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic ; format long ; savepwd = pwd;

% Working directory
cd '..'

mypath = [pwd '/textfiles/']   ;

% Size parameters
run('matlab/parameters.m');

% Load data on price and wage changes
pricehistdatap = importdata([mypath 'data_pdfprices.txt'],' ',0) ;
pricehistdataw = importdata([mypath 'data_pdfwages.txt'],' ',0) ;

% Create plot for data
ip = -0.70:1.40/(2*nump-2):0.70 ; ip = ip' ; xp = [ip';ip']; yp = [pricehistdatap';pricehistdatap'];
iw = -0.70:1.40/(2*numw-2):0.70 ; iw = iw' ; xw = [iw';iw']; yw = [pricehistdataw';pricehistdataw'];

% Load model solution
sol_ss = importdata([mypath '_ss/V10_ss.txt'],' ',0);

% Extract model-generated histograms for versions 2, 3, 4 and 5
ph = sol_ss(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);

distp = ph(1:nump*2-1);
distw = ph(nump*2:nump*2+numw*2-2);

% Create the figure
fig_hist = figure;
subplot(1,2,1)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,distp,'b','LineWidth',1);
    xlabel("Log price changes")
    grid on
    box on
    xlim([-0.5 0.5])
    ylim([0 0.15])
    legend('  Data', '  Model')
    hold off
subplot(1,2,2)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,distw,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off

set(fig_hist,'PaperSize',[20 7],'PaperPosition',[0 0 20 7])
savefig('figures/figs/fig_2.fig')
print(fig_hist,[pwd '/figures/pdfs/fig_2.pdf'],'-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pstep = ip(2)-ip(1);
wstep = iw(2)-iw(1);

ip10 = ip.*0.0;
ip25 = ip.*0.0;
ip05 = ip.*0.0;
for i=1:length(ip)-1
    if ( ip(i) <  -0.025 && ip(i+1) >= -0.025 ) ip25(i) = 1-(abs(ip(i))-0.025)/pstep ; end
    if ( ip(i) >= -0.025 && ip(i)   <=  0.000 ) ip25(i) = 1 ; end
    if ( ip(i) <  -0.100 && ip(i+1) >= -0.100 ) ip10(i) = 1-(abs(ip(i))-0.100)/pstep ; end
    if ( ip(i) >= -0.100 && ip(i)   <=  0.000 ) ip10(i) = 1 ; end
    if ( ip(i) <  -0.050 && ip(i+1) >= -0.050 ) ip05(i) = 1-(abs(ip(i))-0.050)/pstep ; end
    if ( ip(i) >= -0.050 && ip(i)   <=  0.000 ) ip05(i) = 1 ; end
end
massp25_m = sum(ip25.*distp) ;
massp25_d = sum(ip25.*pricehistdatap) ;
massp10_m = sum(ip10.*distp) ;
massp10_d = sum(ip10.*pricehistdatap) ;
massp05_m = sum(ip05.*distp) ;
massp05_d = sum(ip05.*pricehistdatap) ;

iw10 = iw.*0.0;
iw25 = iw.*0.0;
iw05 = iw.*0.0;
for i=1:length(iw)-1
    if ( iw(i) <  -0.025 && iw(i+1) >= -0.025 ) iw25(i) = 1-(abs(iw(i))-0.025)/wstep ; end
    if ( iw(i) >= -0.025 && iw(i)   <=  0.000 ) iw25(i) = 1 ; end
    if ( iw(i) <  -0.100 && iw(i+1) >= -0.100 ) iw10(i) = 1-(abs(iw(i))-0.100)/wstep ; end
    if ( iw(i) >= -0.100 && iw(i)   <=  0.000 ) iw10(i) = 1 ; end
    if ( iw(i) <  -0.050 && iw(i+1) >= -0.050 ) iw05(i) = 1-(abs(iw(i))-0.050)/wstep ; end
    if ( iw(i) >= -0.050 && iw(i)   <=  0.000 ) iw05(i) = 1 ; end
end
massw25_m = sum(iw25.*distw) ;
massw25_d = sum(iw25.*pricehistdataw) ;
massw10_m = sum(iw10.*distw) ;
massw10_d = sum(iw10.*pricehistdataw) ;
massw05_m = sum(iw05.*distw) ;
massw05_d = sum(iw05.*pricehistdataw) ;

fprintf(' \n\n')
fprintf('   ----------------------------------------------  \n')
fprintf('                Changes          Model       Data  \n')
fprintf('   ----------------------------------------------  \n')
fprintf('   Prices   (-0.100, 0.000) %10.4f %10.4f \n', 100*massp10_m,100*massp10_d)
fprintf('            (-0.050, 0.000) %10.4f %10.4f \n', 100*massp05_m,100*massp05_d)
fprintf('            (-0.025, 0.000) %10.4f %10.4f \n', 100*massp25_m,100*massp25_d)
fprintf('   ----------------------------------------------  \n')
fprintf('   Wages    (-0.100, 0.000) %10.4f %10.4f \n', 100*massw10_m,100*massw10_d)
fprintf('            (-0.050, 0.000) %10.4f %10.4f \n', 100*massw05_m,100*massw05_d)
fprintf('            (-0.025, 0.000) %10.4f %10.4f \n', 100*massw25_m,100*massw25_d)
fprintf('   ----------------------------------------------  \n')
fprintf(' \n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load model solutions for versions 2, 3, 4 and 5
sol_ss2 = importdata([mypath '_ss/V15_ss.txt'],' ',0);
sol_ss0 = importdata([mypath '_ss/V12_ss.txt'],' ',0);
sol_ss4 = importdata([mypath '_ss/V13_ss.txt'],' ',0);
sol_ss8 = importdata([mypath '_ss/V14_ss.txt'],' ',0);

% Extract model-generated histograms for versions 2, 3, 4 and 5
ph2 = sol_ss2(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);
ph0 = sol_ss0(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);
ph4 = sol_ss4(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);
ph8 = sol_ss8(4+2*nump*nums+2*numw*numz:3+2*nump*nums+2*numw*numz+2*nump+2*numw);

distp2 = ph2(1:nump*2-1);
distw2 = ph2(nump*2:nump*2+numw*2-2);
distp0 = ph0(1:nump*2-1);
distw0 = ph0(nump*2:nump*2+numw*2-2);
distp4 = ph4(1:nump*2-1);
distw4 = ph4(nump*2:nump*2+numw*2-2);
distp8 = ph8(1:nump*2-1);
distw8 = ph8(nump*2:nump*2+numw*2-2);

ip = -0.70:1.40/(2*nump-2):0.70 ; ip = ip' ; xp = [ip';ip']; yp = [distp';distp'];
iw = -0.70:1.40/(2*numw-2):0.70 ; iw = iw' ; xw = [iw';iw']; yw = [distw';distw'];

% Create the figure
fig_hist = figure;

subplot(4,2,1)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,distp2,'b','LineWidth',1);
    xlabel("Log price changes")
    grid on
    box on
    xlim([-0.5 0.5])
    ylim([0 0.15])
    legend('  2%', ' -2%')
    hold off
subplot(4,2,2)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,distw2,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off

subplot(4,2,3)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,distp0,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    legend('  2%', '  0%')
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off
subplot(4,2,4)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,distw0,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off

subplot(4,2,5)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,distp4,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    legend('  2%', '  4%')
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off
subplot(4,2,6)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,distw4,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off

subplot(4,2,7)
    hold on
    area(xp([2:end end]),yp(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(ip,distp8,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    legend('  2%', '  8%')
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off
subplot(4,2,8)
    hold on
    area(xw([2:end end]),yw(1:end),'FaceColor',[0.8 0.8 0.8],'LineStyle','none','LineWidth',1.0);
    stairs(iw,distw8,'b','LineWidth',1);
    xlabel("Log wage changes")
    grid on
    box on
    xlim([-0.25 0.25])
    ylim([0 0.15])
    hold off

set(fig_hist,'PaperSize',[20 25],'PaperPosition',[0 0 20 25])
savefig('figures/figs/fig_appendix_d1.fig')
print(fig_hist,[pwd '/figures/pdfs/fig_appendix_d1.pdf'],'-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc ; cd(savepwd)


