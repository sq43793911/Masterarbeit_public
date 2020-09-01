%###########################################
% Elementroutine
%###########################################
function [Kte,Me,Be,Qe] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I,It,G,Omega,e)
% Elementroutine: compute Kte, Me

%%
% define empty Kte
Kteux=zeros(12);
Ktev=zeros(12);
Ktevx=zeros(12); 
Ktevxx=zeros(12);
Ktew=zeros(12);
Ktewx=zeros(12); 
Ktewxx=zeros(12);
Ktephi=zeros(12);
Ktephix=zeros(12);
% define empty Me
Me=zeros(12);

Be=zeros(12);
Qe=zeros(12,1);


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
    
   % define N, Nx, Nxx vector linear 
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
        
    % compute Kte, Me, Be, Qe for sampling point of Gauss-integration
    Me=Me + w*( rho*A*( (Nu'*Nu) + (Nv'*Nv) + (Nw'*Nw) ) + rho*I*( (Nvx'*Nvx)+(Nwx'*Nwx)+2*(Nphi'*Nphi) ) )*le/2;
    
    Be=Be+w*( 2*rho*A*Omega*( (Nv'*Nw) - (Nw'*Nv)))*le/2;
    
    Qe=Qe+w*(-rho*A*e*Omega^2*Nv')*le/2;
    
   %% 
%     Kteux=Kteux+w*A*E*(Nux + 2*Nvx*Vx + 2*Nwx*Wx)'*Nux;
%     Ktevx=Ktevx+w*A*E*(Nux*Vx + Nwx*Vx*Wx + Nvx*(Vx^2 + (Ux + Vx^2/2 + Wx^2/2)))'*Nvx;
%     Ktevxx=Ktevxx+w*(E*I*Nvxx)'*Nvxx;
%     Ktewx=Ktewx+w*A*E*(Nux*Wx + Nvx*Vx*Wx + Nwx*(Wx^2 + (Ux + Vx^2/2 + Wx^2/2)))'*Nwx;
%     Ktewxx=Ktewxx+w*(E*I*Nwxx)'*Nwxx;
%     Ktephix=Ktephix+w*(G*It*Nphix)'*Nphix;
    

    %%
    Kteux=Kteux+w*A*E*((1+3*Ux+1.5*Ux^2+0.5*Vx^2+0.5*Wx^2)*Nux+(Vx+Vx*Ux)*Nvx+(Wx+Wx*Ux)*Nwx)'*Nux;
    Ktev=Ktev+w*rho*A*Omega^2*(-Nv'*Nv);
    Ktevx=Ktevx+w*A*E*((Vx+Vx*Ux)*Nux+(Ux+0.5*Ux^2+1.5*Vx^2+0.5*Wx^2)*Nvx+Vx*Wx*Nwx)'*Nvx;
    Ktevxx=Ktevxx+w*E*I*(Nvxx'*Nvxx);
    Ktew=Ktew+w*rho*A*Omega^2*(-Nw'*Nw);
    Ktewx=Ktewx+w*A*E*((Wx+Wx*Ux)*Nux+Vx*Wx*Nvx+(Ux+0.5*Ux^2+0.5*Vx^2+1.5*Wx^2)*Nwx)'*Nwx;
    Ktewxx=Ktewxx+w*E*I*(Nwxx'*Nwxx);
    Ktephi=Ktephi+w*rho*I*2*(-Nphix'*Nphix)*Omega^2;
    Ktephix=Ktephix+w*(G*It*Nphix)'*Nphix;
    
    % matrices Rico
    % nonlinear
    
%     Kteux=Kteux+w*(A*E*(Nux + Nvx*Vx + Nwx*Wx))'*Nux;
%     Ktev=Ktev+w*(-A*Nv*Omega^2*rho)'*Nv;
%     Ktevx=Ktevx+w*((1/2)*A*E*(2*Vx*(Nux + Nwx*Wx) + Nvx*(2*Ux + 3*Vx^2 + Wx^2)))'*Nvx;
%     Ktevxx=Ktevxx+w*(E*I*Nvxx)'*Nvxx;
%     Ktew=Ktew+w*(-A*Nw*Omega^2*rho)'*Nw;
%     Ktewx=Ktewx+w*((1/2)*A*E*(2*(Nux + Nvx*Vx)*Wx + Nwx*(2*Ux + Vx^2 + 3*Wx^2)))'*Nwx;
%     Ktewxx=Ktewxx+w*(E*I*Nwxx)'*Nwxx;
%     Ktephi=Ktephi+w*(-2*I*Nphi*Omega^2*rho)'*Nphi;
%     Ktephix=Ktephix+w*(G*It*Nphix)'*Nphix;
end

Kte=(Kteux+Ktev+Ktevx+Ktevxx+Ktew+Ktewx+Ktewxx+Ktephi+Ktephix)*le/2;
Qe=diag(Qe);
end