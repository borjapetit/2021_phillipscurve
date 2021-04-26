

close all; load GE_dyn_V1_base

shocksize = 1;
xxxcount = 0;

monayshock = [0:0.1:0.10] ;

for i = 1:size(monayshock)

    sumC(i)    = getC(monayshock(i)) ;
    Mshocks(i) = monayshock(i);

    masspchanges = histpchanges(lambda.*ShiftedPdist,Pi);
    masswchanges = histwchanges(rho_MAT.*ShiftedWdist,Pi_w);

    freqsP(xxxcount) = sum(masspchanges);
    freqsW(xxxcount) = sum(masswchanges);

end

figure;

subplot(1,2,1)
plot(100*[0 Mshocks],100*[0 sumC])
xlabel('Size of money shock (%)')
ylabel('Cumulative consumption response (%)')
hold on
axis tight
grid on
subplot(1,2,2)
plot(100*[0 Mshocks],100*[freqpchanges freqsP]);
hold on
plot(100*[0 Mshocks],100*[freqwchanges freqsW]);
axis tight
grid on
legend('Freq. p changes','Freq. w changes')
legend boxoff
xlabel('Size of money shock (%)')
ylabel('Frequency')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outp] = getC(shosksize)

tol      = sqrt(eps) ;
eqcutoff = tol*1000  ;

run([pwd '/matlab/parameters.m']);

sol_ss = importdata([pwd '/textfiles/_ss/V10_ss.txt'],' ',0) ;

run([pwd '/matlab/extract_ss.m']);

load([pwd '/textfiles/_irfs/V10_irf.mat']) ;

rng(1);

TT = 1000;

Shocks      = zeros(nz+nPdist+nWdist+nss,TT);
Shocks(1,:) = shosksize*randn(1,TT);
States      = zeros(nz+nPdist+nWdist+nss,TT);
Jumps       = zeros(nV+nL+nsj,TT);
for time = 2:TT
  States(:,time) = STATEDYNAMICS*States(:,time-1) + Shocks(:,time);
  Jumps(:,time)  = JUMPS*States(:,time);
end
C_path = Cbar + Jumps(nV+nL+1,:);
cons   = (C_path'   - Cbar ) / ( Cbar ) ;

j    = 0;
vec1 = NaN(TT,1);
for i=4:3:size(inflation,1) ; j = j + 1 ;
   vec1(j) = cons(i,1);
end
outp = vec2(2:j) ;
if size(outp,1)>80
    outp = outp(end-79:end);
end
outp = sum(outp)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
