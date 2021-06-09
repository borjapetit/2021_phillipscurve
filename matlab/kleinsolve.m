
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implements Klein's complex QZ decomposition solution method.
%
% Based on Klein's gauss program solve.src, 
%    rewritten for MATLAB 3 March 2008 by Jim Costain.
%    Rewritten using Sylvester equation April 2008.
%    Rewritten to combine solve.src with Klein's MATLAB program solab.m,
%          thus avoiding need to solve Sylvester equation, Jan 2009. 
%
% Idea is to generalize Blanchard/Kahn to cases where their method fails
%    due to invertibility problems. Klein uses QZ decomposition instead
%    of eigendecomposition because it ensures invertibility while decomposing
%    into stable and unstable directions as in Blanchard/Kahn.
%
% Model analyzed is a*E_t(x_t+1) + b*x_t + c*E_t(z_t+1) + d*z_t = 0
%
%    where z_t+1 = phi*z_t + eps_t+1, and eps_t+1 is white noise.
%
%    Here x_t = [k_t ; d_t], where k_t is predetermined: k_t+1 = E_t(k_t+1).
%    User must specify nk, the number of predetermined variables,
%    and must order the x vector accordingly.
%
% Define w_t = [z_t ; k_t].  Klein finds the unique stable solution,
%    if it exists, as follows:
%
%    d_t = F*w_t    w_t+1 = P*w_t + [eps_t+1 ; 0].
%
% Also write entire vector as v_t = [z_t ; k_t ; d_t].
%
% Note: eig<1+eqcutoff is classified as stable.
%   In theory, eqcutoff should be zero.
%   For numerical purposes, may set eqcutoff infinitesimally nonzero.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,P,stableeigs,unstableeigs] = kleinsolve(a,b,c,d,phiMAT,nk,eqcutoff)

nz = size(phiMAT,1);      % # of shocks
nx = size(a,1);           % # of endogenous equations
nd = nx-nk;               % # of endogenous jump variables
nw = nk+nz;               % # of states (endogenous + exogenous)
nv = nx+nz;               % total # of variables

A = [eye(nz) zeros(nz,nx) ;  c  a]; A = full(A);
B = [phiMAT  zeros(nz,nx) ; -d -b]; B = full(B);

save abcd a b c d; clear a b c d

% New representation is:   A * E_t[ z_t+1 ] = B * [ z_t ]
%                                 [ x_t+1 ]       [ x_t ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QZ decomposition:
%   Q*A*Z = S --> Q*A = S*Z'
%   Q*B*Z = T --> Q*B = T*Z'
%   where S and T are upper triangular (or block upper triangular)
%   and Q and Z are unitary, which means Q*Q' = I, Z*Z' = I
%   where prime denotes transpose (if real) or conjugate transpose (if complex).

fprintf('      Start QZ decomposition... ')

[S,T,Q,Z] = qz(A,B);

save AB A B; clear A B;

if ~exist('ordeig')    % For earlier MATLAB versions
   E = eig(T,S);       % Note T,S to compute T_ii/S_ii.    
else
   E = ordeig(T,S);    % Note T,S to compute T_ii/S_ii.
end
selection = abs(E)<=1+eqcutoff;

clear E

tocqz = toc;
fprintf(' %0.1dm %0.1ds \n',floor(tocqz/60),floor(tocqz-60*floor(tocqz/60)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('      Start OrdQZ ... ')

[S,T,Q,Z] = ordqz(S,T,Q,Z,selection); % ~ is Q which isn't used

% This reorders matrices so that |T_ii/S_ii| <= 1+eqcutoff at upper left, and
%                                |T_ii/S_ii|  > 1+eqcutoff at lower right. 

eigenvalues = diag(T)./diag(S);

%  Ratios |T_ii/S_ii| represents generalized eigenvalues of problem:
%     lambda*A*v = B*v  -->  lambda*Q*A*v = Q*B*v  -->  lambda*S*Z'*v = T*Z'*v
%  Define transformed variables y = Z'*v
%     then generalized eigenvalue problem is lambda*S*y = T*y.
%  Note eivals lambda represent rate of increase of system y_t+1 = inv(S)*T*y_t.
%
%  NOTE ORDERING OF EIGENVALUES MEANS --> STABLE BEFORE UNSTABLE 

stableeigs   = eigenvalues(abs(eigenvalues)<=1+eqcutoff);
unstableeigs = eigenvalues(abs(eigenvalues)>1+eqcutoff);
    
tocord = toc-tocqz;

fprintf(' OK: %0.1d seconds \n',floor(tocord));
fprintf('          Number of States: %d    Number of Stable  : %d \n' , nw , length(stableeigs)   ) ;
fprintf('          Number of Jumps : %d    Number of Unstable: %d \n' , nd , length(unstableeigs) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NUMBER OF UNSTABLE EIGENVALUES SHOULD EQUAL nd:
if (length(unstableeigs)<nd || length(stableeigs)<nw)
    fprintf('\n          ERROR!!! Wrong eigenvalue count!!'); fprintf('\n')
    delete abcd.mat
    delete AB.mat
    keyboard 
end

fprintf('      Bulding-up matrices ... ')

Sss = S(1:nw,1:nw);   %stable part: top left block      
Tss = T(1:nw,1:nw);      

Ssu = S(1:nw,nw+1:nv);   
Tsu = T(1:nw,nw+1:nv);   

Suu = S(nw+1:nv,nw+1:nv);   %unstable part: lower right block
Tuu = T(nw+1:nv,nw+1:nv);

% NOW IMPOSE SADDLE PATH STABILITY USING THE QZ DECOMPOSITION.
% System can be rewritten as 
%       Q*A*E_t(v_t+1) = Q*B*v_t --> S*Z'*E_t(v_t+1) = T*Z'*v_t
%       --> E_t(y_t+1) = inv(S)*T*y_t,  where y_t = Z'*v_t
%
% Since S and T are triangular, eigenvalues of inv(S)*T are just T_ii/S_ii.
%
% Break into blocks: y_t = [s_t ; u_t]            -->    STABLE BEFORE UNSTABLE
%    where s_t corresponds to i such that |T_ii/S_ii| <= 1 (stable).
% Then:
%
% [Sss    Ssu] E_t [ s_t+1 ]  =   [Tss  Tsu] * [s_t]
% [ 0     Suu]     [ u_t+1 ]      [ 0   Tuu]   [u_t]
%
% Note dynamics of u_t are independent of s_t, and involve unstable eigenvalues only.
%   So lim j->infinity E_t u_t+j = infinity UNLESS u_t=0.
%
% Thus locally nonexplosive solution requires u_t=0.
%   (NOTE: THIS ORDERING OF STABLE/UNSTABLE IS CRUCIAL FOR THIS SIMPLE SOLUTION.
%          ALWAYS BE CAREFUL ABOUT ORDER WHEN IMPLEMENTING KLEIN METHOD.)
%
% Using u=0, the dynamics become E_t s_t+1 = inv(Sss)*Tss*s_t,
%    which has stable eigenvalues only.

%DYNAMICS of s_t:
sdyn = Sss\Tss;

% NOW TRANSLATE DYNAMICS BACK INTO ORIGINAL VARIABLES.
% By definition, using fact that Z is unitary:
% [w;   =   [z11  z12; * [s;   -->   [s;   =   [z11'  z21'; * [w;
%  d]        z21  z22]    u]          u]        z12'  z22']    d]

% THEREFORE u=0 is equivalent to z12'w + z22'd = 0 -->  d = -inv(z22')*z12'*w
%    Then s = z11'*w + z21'*d = (z11'-inv(z22')*z12'*z21')*w 

% Now using fact that Z is unitary:  
%       Z'*Z=I  -->  z11'*z11+z22'*z22=I  --> z11'+z22'*z22*inv(z11) = inv(z11) 
%       Z'*Z=I  -->  z12'*z11+z22'*z21=0  --> -inv(z22')*z12' = inv(z11)*z21
% 
% can simplify previous formulas:
%       
%       d = -inv(z22')*z12'*w = inv(z11)*z21*w
%       s = (z11'-inv(z22')*z12'*z21')*w = (z11'+inv(z11)*z21*z21')*w = inv(z11)*w
   
z11 = Z(1:nw,1:nw);
z12 = Z(1:nw,nw+1:nv);
z21 = Z(nw+1:nv,1:nw);
z22 = Z(nw+1:nv,nw+1:nv);

tocbuild = toc - tocord - tocqz;
fprintf(' OK: %0.1d seconds \n',floor(tocbuild));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('      Start z11 rank check ... ')

% Getting the rank of z11
[U_z11,SV_z11,V_z11] = svd(z11);   %SVD SHOULD GIVE U*SV*V' = z, with SV diag and U,V unitary.
svdprod1 = U_z11*SV_z11*V_z11';
CHECKSVD1 = max(max(abs(z11 - svdprod1)));

svdprod2 = z11*V_z11;
svdprod3 = U_z11*SV_z11;
CHECKSVD23 = max(max(abs(svdprod2 - svdprod3)));

%SVD is used to calculate rank. So checksvd must be small.
if CHECKSVD1>0.00001 || CHECKSVD23>0.00001
  fprintf('          Problem with check SVD: %d, %d \n',CHECKSVD1,CHECKSVD23)
end

SV_z11 = diag(SV_z11);
tol_rank = max(size(z11)) * eps(max(SV_z11));
rank_z11 = sum(SV_z11 > tol_rank);

if rank_z11  < size(z11,1)
    fprintf('          FATAL ERROR!! Invertibility condition violated \n \n');
    keyboard
else
    fprintf(' OK: '); 
end

tocz11 = toc - tocbuild - tocord - tocqz;
fprintf('%0.1d seconds \n',floor(tocz11));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('      Start consistency checks ... ')

% Relation between d_t and w_t: 
z11i = z11\eye(nw);
F1comp = z21*z11i;
F2comp = -(z22')\z12';
CHECKF = F1comp-F2comp;
checkf = max(max(abs(CHECKF)));

if (checkf >0.00001 )
    fprintf('\n          POSSIBLE ERROR: The two formulas for F give different answers! %d \n          ',checkf)
end

% Dynamics of w_t: 
%   E_t w_t+1 = z11 E_t s_t+1 = z11*inv(Sss)*Tss*s_t = z11*inv(Sss)*Tss*inv(z11)*w_t
Pcomp = z11*sdyn*z11i;

maximagF = max(max(abs(imag(F1comp))));
if maximagF>0.00001
    fprintf('\n          Warning! Large imaginary part observed in F: %d \n         ', maximagF)
end
maximagP = max(max(abs(imag(Pcomp))));
if maximagP>0.00001
    fprintf('\n          Warning! Large imaginary part observed in P: %d \n         ',maximagP)
end
biggeststable = max(abs(eig(Pcomp)));
if biggeststable > 1+eqcutoff
    fprintf('\n          FATAL ERROR!! large eigenvector classified as stable!! \n \n'); fprintf('\n')
    keyboard
elseif biggeststable > 1
   Pcomp=Pcomp/biggeststable;  % Note: biggeststable infinitesimally greater than one: rescaling down
end

F = real(F1comp);
P = real(Pcomp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTATION IS FINISHED.  Now check identities implied by the model.

Pz = P(1:nz,:);
Pk = P(nz+1:end,:);

% KLEIN SETUP IS: a*E_t(x_t+1) + b*x_t + c*E_t(z_t+1) + d*z_t = 0
% starting from any arbitrary w_t=[z_t;k_t], which determines d_t, etc.
%
% Here E_t x_t+1 = [k_t+1; = [ Pk; * [z_t;  and x_t = [k_t; = [0 I; * [z_t; 
%                   d_t+1]    F*P]    k_t]             d_t]     F ]    k_t]
%
% Similarly, A*E_t(v_t+1) = B*v_t from arbitrary w_t=[z_t;k_t], but not
% from arbitrary v_t.
%
% IF SOLUTION IS CORRECT THE FOLLOWING MATRICES SHOULD BE IDENTICALLY ZERO:

load abcd; CHECKKLEINSOLVE1 = a*[Pk ; F*P]  + b*[[zeros(nk,nz) eye(nk)] ; F] + c*Pz + d*[eye(nz) zeros(nz,nk)]; clear a b c d;
load AB;   CHECKKLEINSOLVE2 = A*[P ; F*P]   - B*[eye(nw) ; F]; clear A B;
           CHECKKLEINSOLVE3 = S*Z'*[P; F*P] - T*Z'*[eye(nw) ; F];

checkklein1 = max(max(abs(real(CHECKKLEINSOLVE1))));
checkklein2 = max(max(abs(real(CHECKKLEINSOLVE2))));
checkklein3 = max(max(abs(real(CHECKKLEINSOLVE3))));

if (checkklein1 > 0.00001 || checkklein2 > 0.00001 || checkklein3 > 0.00001 )
    fprintf(' POSSIBLE ERROR?\n');
    fprintf('         Identities of Klein model not satisfied?\n');
    fprintf('         %d = %d = %d ??? \n',checkklein1,checkklein2,checkklein3);
else 
    fprintf(' OK: ');
    toccon = toc - tocz11 - tocbuild - tocord - tocqz;
    fprintf(' %0.1d seconds \n',floor(toccon));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear toccon tocz11 tocord tocqz tocbuild

delete abcd.mat
delete AB.mat
