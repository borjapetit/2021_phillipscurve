
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE GENERATES THE FIGURE 8 IN THE PAPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; close all ; clc ; tic ; format long ; savepwd = pwd;

% Working directory
cd '..'

% Add paths
addpath([pwd '/matlab/'])
addpath([pwd '/figures/'])

% Load model parameters
sol_ss = importdata([pwd '/textfiles/_ss/V10_ss.txt'],' ',0);
run('matlab/parameters.m');
run('matlab/extract_ss.m');

clearvars -except nums nump numz numw savepwd mu

% Load price and wage setting choices
lambda = importdata([pwd '/textfiles/_ss/_V10_lambda_ss.txt'],' ',0); lambda = reshape(lambda   ,nums,nump) ; lambda = lambda' ;
Pi     = importdata([pwd '/textfiles/_ss/_V10_pi_ss.txt'    ],' ',0); Pi     = reshape(Pi       ,nums,nump) ; Pi     = Pi'     ;
rho    = importdata([pwd '/textfiles/_ss/_V10_rho_ss.txt'   ],' ',0); rho    = reshape(rho      ,numz,numw) ; rho    = rho'    ;
Pi_w   = importdata([pwd '/textfiles/_ss/_V10_piw_ss.txt'   ],' ',0); Pi_w   = reshape(Pi_w,numw,numz,numw) ; Pi_w   = permute(Pi_w,[3,2,1]) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0 ; shocksize = 1;

fprintf('\n') 
fprintf('\n') 

for mstep = 0:0.4:12; i = i + 1;
    
    fprintf(' Money shock = %4.2f%% \n',mstep)
    
    [CC(i),MM(i),SPdist,SWdist] = getirf(mstep);
    freqsP(i) = sum(sum(lambda.*SPdist)) ;
    freqsW(i) = sum(sum(rho.*SWdist))    ;

end

f1 = figure;
subplot(1,2,1)
    hold on
    box on
    plot(100*MM,100*CC,'-o','LineWidth',1.0,'MarkerSize',5)
    plot(100*MM,100*CC.*0,'-black','LineWidth',1.0,'MarkerSize',5)
    ylim([-1 3.5])
    xlim([100*MM(1) 100*MM(end)])
    xlabel('Size of money shock (%)')
    ylabel('Cum. consumption response (%)')
    grid on
subplot(1,2,2)
    hold on
    box on
    plot(100*MM,100*freqsP,'-d','LineWidth',1.0,'MarkerSize',5)
    plot(100*MM,100*freqsW,'-s','LineWidth',1.0,'MarkerSize',5)
    ylim([0 100])
    xlim([100*MM(1) 100*MM(end)])
    grid on
    legend('Freq. price changes','Freq. wage changes','Location','NorthWest')
    legend boxoff
    xlabel('Size of money shock (%)')
    ylabel('Frequency')

set(f1,'PaperSize',[20 7],'PaperPosition',[0 0 20 7])
savefig('figures/figs/fig_8.fig')
print(f1,'figures/pdfs/fig_8','-dpdf')

fprintf('\n') ; toc ; fprintf('\n')






cd(savepwd)

function [CC,MM,SPdist,SWdist] = getirf(sstep)

    sol_ss = importdata([pwd '/textfiles/_ss/V10_ss.txt'],' ',0);
    
    run('matlab/parameters.m');
    run('matlab/extract_ss.m');

    load([pwd '/textfiles/_irfs/V10_irf.mat']) ;

    clear STATEHISTORY ; clear JumpHistory ; TT = 24 ;

    Shocks   = zeros(nz+nPdist+nWdist+nss,TT) ;
    States   = zeros(nz+nPdist+nWdist+nss,TT) ;
    Jumps    = zeros(nV+nL+nsj,TT)            ;
    Inf_path = zeros(1,TT)                    ;
    C_path   = zeros(1,TT)                    ;
    prices   = zeros(1,TT)                    ;
    M_path   = zeros(1,TT)                    ;
    
    SPdist = Rmatrix(nump,exp(sstep/100)+eps,pstep)*Pdist;
    SWdist = Rmatrix(numw,exp(sstep/100)+eps,wstep)*Wdist;
    
    Pvector = SPdist-Pdist ; Pvector = Pvector';%zeros(size(Pdist,1)*size(Pdist,2));
    Wvector = SWdist-Wdist ; Wvector = Wvector';%zeros(size(Wdist,1)*size(Wdist,2));

    States(:,1) = [zeros(nz,1);Pvector(:);Wvector(:);zeros(nss,1)];
    for time = 2:TT
      States(:,time) = STATEDYNAMICS*States(:,time-1) + Shocks(:,time);
    end
    for time = 1:TT
      Jumps(:,time) = JUMPS*States(:,time);
    end
    
    Inf_path    = mu          + States(end,:);
    Inf_path(1) = Inf_path(1) + log( exp(sstep/100) + eps );
 
    C_path = Cbar + Jumps(nV+nL+1,:) ;
    prices = cumsum(Inf_path - 1);
    M_path = ones(1,TT)*log(mu) ;
    
    CC = sum( C_path-Cbar ) ;
    MM = prices(end)-sum(M_path(2:end)) - log(mu) ;
    
end



function R = Rmatrix(numpp,infl,pstepp)
  R = sparse(numpp,numpp); 
  nowoffset = log(infl)/pstepp;
  if nowoffset==0
      R = eye(numpp);
  else
      remoffset = nowoffset - floor(nowoffset);
      R(1,1:ceil(nowoffset)) = 1;             
      startfirstdiag = [max([1 ; -floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];         
      endfirstdiag   = [min([numpp-1 ; numpp-ceil(nowoffset)])  min([numpp ; numpp+floor(nowoffset)])]; 
      startsecdiag = [max([2 ; 1-floor(nowoffset)])  max([1 ; 1+ceil(nowoffset)])];               
      endsecdiag   = [min([numpp ; numpp-ceil(nowoffset)+1])  min([numpp ; numpp+floor(nowoffset)])];   
      R(startfirstdiag(1):endfirstdiag(1),startfirstdiag(2):endfirstdiag(2)) = ...
           remoffset*speye(numpp-ceil(abs(nowoffset)));
      R(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) = ...
      R(startsecdiag(1):endsecdiag(1),startsecdiag(2):endsecdiag(2)) + ...
           (1-remoffset)*speye(numpp-ceil(abs(nowoffset)));
      R(numpp,numpp+floor(nowoffset)+1:numpp) = 1; 
  end
end

function masspchangesvec=histpchanges(Adj,PMat)
    run('matlab/parameters.m');
    masspchanges = zeros(2*nump-1,nums);
    for i = 1:nums
    adjprobMAT = Adj(:,i)*(PMat(:,i))';
      for j = 1-nump:1:nump-1
        masspchanges(j+nump,i) = sum(diag(adjprobMAT,j));
      end
    end
    masspchangesvec = sum(masspchanges,2);
end

function masspchangesvec=histwchanges(Adj,PMat)
    run('matlab/parameters.m');
    masspchanges = zeros(2*numw-1,numz);
    adjprobMAT   = zeros(numw,numw);
    for i = 1:numz
        for k = 1:nump
            adjprobMAT(k,:) = Adj(k,i).*PMat(k,i,:);
        end
        for j = 1-numw:1:numw-1
            masspchanges(j+numw,i) = sum(diag(adjprobMAT,j));
        end
    end
    masspchangesvec = sum(masspchanges,2);
end






