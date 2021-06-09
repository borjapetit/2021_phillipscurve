
sol_dyn = sol_dyn';

ix      = 1:Ntot;
ixnext  = Ntot+1:2*Ntot;
iz      = 2*Ntot+1:2*Ntot+nz;
iznext  = 2*Ntot+nz+1:2*(Ntot+nz);

ajac    = sol_dyn(:,ixnext) ;
bjac    = sol_dyn(:,ix);
cjac    = sol_dyn(:,iznext) ;
djac    = sol_dyn(:,iz);

clear ix ixnext iz iznext