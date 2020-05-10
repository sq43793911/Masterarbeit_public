%###########################################
% Elementroutine
%###########################################
function [Kte,Me,Qe] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I,q)
% Elementroutine: compute Kte, Me

%%
% % define empty Kte
Kteux=zeros(10,10); 
Ktevx=zeros(10,10); 
Ktevxx=zeros(10,10); 
Ktewx=zeros(10,10); 
Ktewxx=zeros(10,10); 
% % define empty Me
Me=zeros(10,10);
Q=zeros(4,1);
Qe=zeros(10,1);

Nq1=zeros(4,1);
Nq2=zeros(2,1);


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
    N=[ 0.5-xi/2, 0,                     0,                          0,                     0,                          0.5+xi/2, 0,                       0,                          0,                      0; ...
        0,        1/2-(3*xi)/4+(xi^3)/4, 1/4-xi/4-(xi^2)/4+(xi^3)/4, 0,                     0,                          0,        1/2+(3*xi)/4-(xi^3)/4,  -1/4-xi/4+(xi^2)/4+(xi^3)/4, 0,                      0; ...
        0,        -3/4+(3*xi^2)/4,       -1/4-xi/2+(3*xi^2)/4,       0,                     0,                          0,        3/4-(3*xi^2)/4,         -1/4+xi/2+(3*xi^2)/4,        0,                      0; ...
        0,        0,                     0,                          1/2-(3*xi)/4+(xi^3)/4, 1/4-xi/4-(xi^2)/4+(xi^3)/4, 0,        0,                      0,                           1/2+(3*xi)/4-(xi^3)/4, -1/4-xi/4+(xi^2)/4+(xi^3)/4;
        0,        0,                     0,                         -3/4+(3*xi^2)/4,       -1/4-xi/2+(3*xi^2)/4,        0,        0,                      0,                           3/4-(3*xi^2)/4,        -1/4+xi/2+(3*xi^2)/4];
    
    Nx=[ -0.5, 0,               0,                    0,                0,                    0.5,  0,                0,                    0,                 0; ...
          0,  -3/4+(3*xi^2)/4, -1/4-xi/2+(3*xi^2)/4,  0,                0,                    0,    3/4-(3*xi^2)/4,  -1/4+xi/2+(3*xi^2)/4,  0,                 0; ...
          0,   3*xi/2,         -1/2+3*xi/2,           0,                0,                    0,   -3*xi/2,           1/2+3*xi/2,           0,                 0; ...
          0,   0,               0,                   -3/4+(3*xi^2)/4,  -1/4-xi/2+(3*xi^2)/4,  0,    0,                0,                    3/4-(3*xi^2)/4,   -1/4+xi/2+(3*xi^2)/4; ...
          0,   0,               0,                    3*xi/2,           -1/2+3*xi/2,          0,    0,                0,                    -3*xi/2,           1/2+3*xi/2]*(2/le);
    
    Nxx=[ 0,  0,       0,          0,       0,          0,  0,       0,           0,      0; ...
          0,  3*xi/2, -1/2+3*xi/2, 0,       0,          0, -3*xi/2,  1/2+3*xi/2,  0,      0; ...
          0,  3/2,     3/2,        0,       0,          0, -3/2,     3/2,         0,      0; ...
          0,  0,       0,          3*xi/2, -1/2+3*xi/2, 0,  0,       0,          -3*xi/2, 1/2+3*xi/2; ...
          0,  0,       0,          3/2,     3/2,        0,  0,       0,          -3/2,    3/2]*((2/le)^2);  
      
    Nu=N(1,:);
    Nw=N(2,:);
    Nv=N(4,:);
    %%
    Nq1(1)=Nw(2);
    Nq1(2)=Nw(3);
    Nq1(3)=Nw(7);
    Nq1(4)=Nw(8);
    
    Nq2(1)=Nu(1);
    Nq2(2)=Nu(6);
    
    Nq=Nq1*Nq2';
    
    %%
    Nux=Nx(1,:);
    Nwx=Nx(2,:);
    Nvx=Nx(4,:);
    
    Nuxx=Nxx(1,:);
    Nwxx=Nxx(2,:);
    Nvxx=Nxx(4,:);
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me +  rho * A * (Nu' * Nu) * w*le/2 + rho * A * (Nw' * Nw) * w*le/2 + rho * A * (Nv' * Nv) * w *le/2 ;
    
      
    Kteux=Kteux+w*(A*E*Nux + 2*A*E*Nvx*Vx + 2*A*E*Nwx*Wx)'*Nux*le/2;
    Ktevx=Ktevx+w*(A*E*Nux*Vx + A*E*Nwx*Vx*Wx + Nvx*(A*E*Vx^2 + A*E*(Ux + Vx^2/2 + Wx^2/2)))'*Nvx*le/2;
    Ktevxx=Ktevxx+w*(E*I*Nvxx)'*Nvxx*le/2;
    Ktewx=Ktewx+w*(A*E*Nux*Wx + A*E*Nvx*Vx*Wx + Nwx*(A*E*Wx^2 + A*E*(Ux + Vx^2/2 + Wx^2/2)))'*Nwx*le/2;
    Ktewxx=Ktewxx+w*(E*I*Nwxx)'*Nwxx*le/2;
    
    Q=Q+le*Nq*[q;q]*w;
end
Kte=Kteux+Ktevx+Ktevxx+Ktewx+Ktewxx;

Qe(1)=0;
Qe(2)=Q(1);
Qe(3)=Q(2);
Qe(4)=0;
Qe(5)=0;

Qe(6:end)=Qe(1:5);

Qe=diag(Qe);

end