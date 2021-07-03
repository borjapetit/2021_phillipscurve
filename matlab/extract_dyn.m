
%sol_dyn = sol_dyn';
sol_dyn  = reshape(sol_dyn,[9047,18100]);

ix      = 1:Ntot;
ixnext  = Ntot+1:2*Ntot;
iz      = 2*Ntot+1:2*Ntot+nz;
iznext  = 2*Ntot+nz+1:2*(Ntot+nz);

ajac    = sol_dyn(:,ixnext) ;
bjac    = sol_dyn(:,ix) ;
cjac    = sol_dyn(:,iznext) ;
djac    = sol_dyn(:,iz) ;

fprintf(' Check matrix size... \n');
fprintf('    Size solve_dyn: %1.0f x %1.0f \n',size(sol_dyn,1),size(sol_dyn,2))
fprintf('    Size ajac:      %1.0f x %1.0f \n',size(ajac,1),size(ajac,2))
fprintf('    Size bjac:      %1.0f x %1.0f \n',size(bjac,1),size(bjac,2))
fprintf('    Size cjac:      %1.0f x %1.0f \n',size(cjac,1),size(cjac,2))
fprintf('    Size djac:      %1.0f x %1.0f \n',size(djac,1),size(djac,2))

clear ix ixnext iz iznext