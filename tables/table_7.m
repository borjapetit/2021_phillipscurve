

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 5 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic; format long

cd '/Users/borjapetit/Dropbox/projects/2016_lpw/code_public/codes/'

mypath = [pwd '/textfiles/_irfs/']   ;

PC = ones(4,5);
CL = ones(4,5);

eval(['load ' pwd '/textfiles/_irfs/V15_irf.mat irf' ] ); PC(1,1) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V35_irf.mat irf' ] ); PC(2,1) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V55_irf.mat irf' ] ); PC(3,1) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V65_irf.mat irf' ] ); PC(4,1) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V12_irf.mat irf' ] ); PC(1,3) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V32_irf.mat irf' ] ); PC(2,3) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V52_irf.mat irf' ] ); PC(3,3) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V62_irf.mat irf' ] ); PC(4,3) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V10_irf.mat irf' ] ); PC(1,5) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V30_irf.mat irf' ] ); PC(2,5) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V50_irf.mat irf' ] ); PC(3,5) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V60_irf.mat irf' ] ); PC(4,5) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V13_irf.mat irf' ] ); PC(1,6) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V33_irf.mat irf' ] ); PC(2,6) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V53_irf.mat irf' ] ); PC(3,6) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V63_irf.mat irf' ] ); PC(4,6) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V14_irf.mat irf' ] ); PC(1,7) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V34_irf.mat irf' ] ); PC(2,7) = sum(irf.pinflation)/sum(irf.labor); 
eval(['load ' pwd '/textfiles/_irfs/V54_irf.mat irf' ] ); PC(3,7) = sum(irf.pinflation)/sum(irf.labor);
eval(['load ' pwd '/textfiles/_irfs/V64_irf.mat irf' ] ); PC(4,7) = sum(irf.pinflation)/sum(irf.labor); 

fprintf('\n')
fprintf('\n')
fprintf('    TABLE 7. PHILLIPS MULTIPLIERS AT DIFFERENT INFLATION RATES\n')
fprintf('             AND NOISE PARAMETERS')
fprintf('\n')
fprintf('\n')
fprintf('    **************************************************** \n')
fprintf('                         Base       FP       FW     FPFW \n')
fprintf('    **************************************************** \n')
fprintf('    Inflation = -2%%    %6.3f   %6.3f   %6.3f   %6.3f    \n', PC(:,1))
fprintf('    Inflation =  0%%    %6.3f   %6.3f   %6.3f   %6.3f    \n', PC(:,3))
fprintf('    Inflation =  2%%    %6.3f   %6.3f   %6.3f   %6.3f    \n', PC(:,5))
fprintf('    Inflation =  4%%    %6.3f   %6.3f   %6.3f   %6.3f    \n', PC(:,6))
fprintf('    Inflation =  8%%    %6.3f   %6.3f   %6.3f   %6.3f    \n', PC(:,7))
fprintf('    **************************************************** \n')
fprintf('\n\n\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PClr  = zeros(5,1) ; 
PCsrp = zeros(5,1) ; 
PCsrw = zeros(5,1) ; 

j = 0;
for version_kappas=[15 12 10 13 14] ; j = j + 1 ;
    sol_irf = importdata([pwd '/textfiles/_irfs/V' num2str(version_kappas) '_irf.mat'],' ',0);   
    slope_pc(j,1)  = sum(sol_irf.pinflation)/sum(sol_irf.labor) ;
end

xx = [1 2 3 4 5] ;

f = figure;
    hold on;
    bar(xx,slope_pc')
    xticks([1 2 3 4 5])
    xticklabels({'-2%','0%','2%','4%','8%'})
    ylim([0 0.5])
    xlabel('Inflation rate')
    title('Phillips Multiplier')
    hold off

set(f,'PaperSize',[10 7],'PaperPosition',[0 0 10 7])
savefig('figures/_others/fig_multiplier.fig')
print(f,[pwd '/figures/_others/fig_multiplier.pdf'],'-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




