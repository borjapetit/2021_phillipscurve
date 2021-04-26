

pPdist = 3+1:nPdist+3;
pV     = nPdist+3+1:nPdist+3+nV;
pWdist = nPdist+3+nV+1:nPdist+3+nV+nWdist;
pL     = nPdist+3+nV+nWdist+1:nPdist+3+nV+nWdist+nL;

Cbar   = sol_ss(1);
Nbar   = sol_ss(2);
wbar   = sol_ss(3);
Pdist  = sol_ss(pPdist); 
V      = sol_ss(pV); 
Wdist  = sol_ss(pWdist);
L      = sol_ss(pL); 
histmo = sol_ss(4+2*nPdist+2*nWdist:3+2*nPdist+2*nWdist+2*nump+2*numw);
istat  = sol_ss(3+2*nPdist+2*nWdist+2*nump+2*numw+1:end);

mbar    = nu*Cbar^gamma/(1-1/Rss);
Pdist   = reshape(Pdist   ,nums,nump) ; Pdist  = Pdist'  ;
V       = reshape(V       ,nums,nump) ; V      = V'      ;
Wdist   = reshape(Wdist   ,numz,numw) ; Wdist  = Wdist'  ;
L       = reshape(L       ,numz,numw) ; L      = L'      ;

hist_p = histmo(1:nump*2-1) ; 
hist_w = histmo(nump*2:nump*2+numw*2-2) ;

istat  = sol_ss(3+2*nPdist+2*nWdist+2*nump+2*numw+1:end);

stat.freq_p     = istat(1)  ; stat.freq_w     = istat(1+length(istat)/2) ;
stat.avech_p    = istat(2)  ; stat.avech_w    = istat(2+length(istat)/2) ;
stat.aveabsch_p = istat(3)  ; stat.aveabsch_w = istat(3+length(istat)/2) ;
stat.std_p      = istat(4)  ; stat.std_w      = istat(4+length(istat)/2) ;
stat.skech_p    = istat(5)  ; stat.skech_w    = istat(5+length(istat)/2) ;
stat.kurch_p    = istat(6)  ; stat.kurch_w    = istat(6+length(istat)/2) ;
stat.ch0_p      = istat(7)  ; stat.ch0_w      = istat(7+length(istat)/2) ;
stat.absch5_p   = istat(8)  ; stat.absch5_w   = istat(8+length(istat)/2) ;
stat.absch25_p  = istat(9)  ; stat.absch25_w  = istat(9+length(istat)/2) ;
stat.losses_p   = istat(10) ; stat.losses_w   = istat(10+length(istat)/2) ;
stat.errors_p   = istat(11) ; stat.errors_w   = istat(11+length(istat)/2) ;
stat.setting_p  = istat(12) ; stat.setting_w  = istat(12+length(istat)/2) ;
stat.timing_p   = istat(13) ; stat.timing_w   = istat(13+length(istat)/2) ;
stat.totalloss  = istat(14) ;

clear pricehistmodel pPdist pV pWdist pL istat

