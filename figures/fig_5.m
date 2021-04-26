
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 5 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

% Working directory
cd '/Users/borjapetit/Dropbox/projects/2016_lpw/code_public/codes/'

% Print figures for versions "base" "FP", "FW" and "FPFW"
irfsall = figure; set(irfsall ,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
eval(['load ' pwd '/textfiles/_irfs/V10_irf.mat irf']) ; linescolor = '-o' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V30_irf.mat irf']) ; linescolor = '-^' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V50_irf.mat irf']) ; linescolor = '-d' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V60_irf.mat irf']) ; linescolor = '-s' ; run('matlab/plot_irf.m')
subplot(2,3,1); legend(' Baseline' ,' (FP) Flexi. prices' ,' (FW) Flexi. wages',' (FPFW) Both flexible')
savefig('figures/figs/fig_5.fig')
print(irfsall,'figures/pdfs/fig_5','-dpdf')
clear irfsall
