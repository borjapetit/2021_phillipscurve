
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES TABLE 6 IN THE PAPER
% It finds the best fiiting Calvo Parametrs for 4 different inflation rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic ; format long ; savepwd = pwd;

% Working directory
cd '..'

% Path to Dynare
addpath('/Applications/Dynare/4.6.4/matlab/')

path0 = pwd;

fprintf('\n')

% Load solution for inflation rate of 2%
clearvars -except path0 savepwd
eval( ['load ' pwd '/textfiles/_irfs/V10_irf' ] );
eval( ['save ' pwd '/matlab/V10_irf'  ] );

% Load solution for inflation rate of 0%
clearvars -except path0 savepwd
eval( ['load ' pwd '/textfiles/_irfs/V12_irf' ] );
eval( ['save ' pwd '/matlab/V12_irf'  ] );

% Load solution for inflation rate of 4%
clearvars -except path0 savepwd
eval( ['load ' pwd '/textfiles/_irfs/V13_irf' ] );
eval( ['save ' pwd '/matlab/V13_irf'  ] );

% Load solution for inflation rate of 8%
clearvars -except path0 savepwd
eval( ['load ' pwd '/textfiles/_irfs/V14_irf' ] );
eval( ['save ' pwd '/matlab/V14_irf'  ] );

% Find calvo parameters
cd([pwd '/matlab/'])
dynare calvow
calvopars = 1-x_opt_hat';

% Clear workspace
clearvars -except path0 savepwd calvopars
delete 'V10_irf.mat' 'V12_irf.mat' 'V13_irf.mat' 'V14_irf.mat' 'calvow_results.mat' 'calvow.log' 'H.dat' 'g1.mat' ;
rmdir('calvow', 's')
rmdir('+calvow', 's')

% Print table
clc
cd '..'
diary 'tables/table_6.txt'
fprintf('\n')
fprintf('\n')
fprintf('    TABLE 6. BEST FITTING CALVO MODELS')
fprintf('\n')
fprintf('\n')
fprintf('    ********************************** \n')
fprintf('                       PRICES    WAGES \n')
fprintf('    ********************************** \n')
fprintf('    Inflation =  0%%    %6.4f   %6.4f  \n', calvopars(2,:))
fprintf('    Inflation =  2%%    %6.4f   %6.4f  \n', calvopars(1,:))
fprintf('    Inflation =  4%%    %6.4f   %6.4f  \n', calvopars(3,:))
fprintf('    Inflation =  8%%    %6.4f   %6.4f  \n', calvopars(4,:))
fprintf('    ********************************** \n')
fprintf('\n\n\n\n')
diary off

toc ; cd(savepwd)



