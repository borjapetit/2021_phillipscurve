
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 8 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

% Working directory
cd ' %% Your Own Directory '

% Print figure - paper version (V1, V3, V5 and V6)
irfsall = figure; set(irfsall ,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
eval(['load ' pwd '/textfiles/_irfs/V15_irf.mat irf']) ; linescolor = '-o' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V12_irf.mat irf']) ; linescolor = '-s' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V10_irf.mat irf']) ; linescolor = '-d' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V13_irf.mat irf']) ; linescolor = '-^' ; run('matlab/plot_irf.m')
eval(['load ' pwd '/textfiles/_irfs/V14_irf.mat irf']) ; linescolor = '-v' ; run('matlab/plot_irf.m')
subplot(2,3,1); legend(' -2%','  0%' ,'  2%','  4%','  8%')
savefig('figures/figs/fig_6.fig')
print(irfsall,'figures/pdfs/fig_6','-dpdf')
clear irfsall
