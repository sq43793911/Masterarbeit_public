%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I,It,G)
% Elementroutine: compute Kte, Me

%%
% % define empty Kte
Kteux=zeros(12); 
Ktevx=zeros(12); 
Ktevxx=zeros(12); 
Ktewx=zeros(12); 
Ktewxx=zeros(12);
Ktephix=zeros(12);
% % define empty Me
Me=zeros(12);


%%
% xiVec=[-sqrt(1/3),sqrt(1/3)];   % define sampling points for Gauss-quadrature      
% wVec =[1,1];   % weights for sampling points of Gauss-quadrature 

%
% xiVec=[-sqrt(3/5),0,sqrt(3/5)];  % define sampling points for Gauss-quadrature
% wVec=[5/9,8/9,5/9];             % weights for sampling points of Gauss-quadrature

xiVec=[-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
wVec=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
%%

for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
   % define N, B vector linear 
    N=[ 0.5-xi/2, 0,                     0,                          0,                     0,                          0,        0.5+xi/2, 0,                       0,                          0,                      0,                          0; ...
        0,        1/2-(3*xi)/4+(xi^3)/4, 1/4-xi/4-(xi^2)/4+(xi^3)/4, 0,                     0,                          0,        0,        1/2+(3*xi)/4-(xi^3)/4,  -1/4-xi/4+(xi^2)/4+(xi^3)/4, 0,                      0,                          0; ...
        0,        -3/4+(3*xi^2)/4,       -1/4-xi/2+(3*xi^2)/4,       0,                     0,                          0,        0,        3/4-(3*xi^2)/4,         -1/4+xi/2+(3*xi^2)/4,        0,                      0,                          0; ...
        0,        0,                     0,                          1/2-(3*xi)/4+(xi^3)/4, 1/4-xi/4-(xi^2)/4+(xi^3)/4, 0,        0,        0,                      0,                           1/2+(3*xi)/4-(xi^3)/4, -1/4-xi/4+(xi^2)/4+(xi^3)/4, 0; ...
        0,        0,                     0,                         -3/4+(3*xi^2)/4,       -1/4-xi/2+(3*xi^2)/4,        0,        0,        0,                      0,                           3/4-(3*xi^2)/4,        -1/4+xi/2+(3*xi^2)/4,        0; ...
        0,        0,                     0,                          0,                     0,                          0.5-xi/2, 0,        0,                      0,                           0,                      0,                          0.5+xi/2];
    
    Nx=[ -0.5, 0,               0,                    0,                0,                    0,   0.5,  0,                0,                    0,                 0,                   0; ...
          0,  -3/4+(3*xi^2)/4, -1/4-xi/2+(3*xi^2)/4,  0,                0,                    0,   0,    3/4-(3*xi^2)/4,  -1/4+xi/2+(3*xi^2)/4,  0,                 0,                   0; ...
          0,   3*xi/2,         -1/2+3*xi/2,           0,                0,                    0,   0,   -3*xi/2,           1/2+3*xi/2,           0,                 0,                   0; ...
          0,   0,               0,                   -3/4+(3*xi^2)/4,  -1/4-xi/2+(3*xi^2)/4,  0,   0,    0,                0,                    3/4-(3*xi^2)/4,   -1/4+xi/2+(3*xi^2)/4, 0; ...
          0,   0,               0,                    3*xi/2,          -1/2+3*xi/2,           0,   0,    0,                0,                   -3*xi/2,            1/2+3*xi/2           0; ...
          0,   0,               0,                    0,                0,                   -0.5, 0,    0,                0,                    0,                 0,                   0.5]*(2/le);
    
    Nxx=[ 0,  0,       0,          0,       0,          0, 0,  0,       0,           0,      0,               0; ...
          0,  3*xi/2, -1/2+3*xi/2, 0,       0,          0, 0, -3*xi/2,  1/2+3*xi/2,  0,      0,               0; ...
          0,  3/2,     3/2,        0,       0,          0, 0, -3/2,     3/2,         0,      0,               0; ...
          0,  0,       0,          3*xi/2, -1/2+3*xi/2, 0, 0,  0,       0,          -3*xi/2, 1/2+3*xi/2,      0; ...
          0,  0,       0,          3/2,     3/2,        0, 0,  0,       0,          -3/2,    3/2,             0; ...
          0,  0,       0,          0,       0,          0, 0,  0,       0,           0,      0,               0]*((2/le)^2);
      
    Nu=N(1,:);
    Nw=N(2,:);
    Nv=N(4,:);
    Nphi=N(6,:);
    
    Nux=Nx(1,:);
    Nwx=Nx(2,:);
    Nvx=Nx(4,:);
    Nphix=Nx(6,:);
    
    Nuxx=Nxx(1,:);
    Nwxx=Nxx(2,:);
    Nvxx=Nxx(4,:);
    Nphixx=Nxx(6,:);
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me +  rho * A * (Nu' * Nu) * w + rho * A * (Nw' * Nw) * w + rho * A * (Nv' * Nv) * w + rho * A * (Nphi' * Nphi) * w;
    
      
    Kteux=Kteux+w*(A*E*Nux + 2*A*E*Nvx*Vx + 2*A*E*Nwx*Wx)'*Nux;
    Ktevx=Ktevx+w*(A*E*Nux*Vx + A*E*Nwx*Vx*Wx + Nvx*(A*E*Vx^2 + A*E*(Ux + Vx^2/2 + Wx^2/2)))'*Nvx;
    Ktevxx=Ktevxx+w*(E*I*Nvxx)'*Nvxx;
    Ktewx=Ktewx+w*(A*E*Nux*Wx + A*E*Nvx*Vx*Wx + Nwx*(A*E*Wx^2 + A*E*(Ux + Vx^2/2 + Wx^2/2)))'*Nwx;
    Ktewxx=Ktewxx+w*(E*I*Nwxx)'*Nwxx;
    Ktephix=Ktephix+w*(G*It*Nphix)'*Nphix;
end
Kte=(Kteux+Ktevx+Ktevxx+Ktewx+Ktewxx+Ktephix)*le/2;
end