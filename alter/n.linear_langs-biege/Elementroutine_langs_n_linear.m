%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_langs_n_linear(A,E,rho,le,Ux,Vx,Wx,I)
% Elementroutine: compute Kte, Me

%%
% % define empty Kte
Kte=zeros(10,10); 
% % define empty Me
Me=zeros(10,10);

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
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + rho * A * (N' * N) * w*le/2;
    Kte=Kte + E * A * (2+6*Ux+3*Ux^2+Vx^2+Wx^2 + 2*Vx+2*Vx*Ux + 2*Wx+2*Wx*Ux ...
        + 2*Vx+2*Vx*Ux + 2*Ux+Ux^2+3*Vx^2+Wx^2 + 2*Wx+2*Vx*Wx ...
        + 2*Wx+2*Wx*Ux + 2*Vx*Wx + 2*Ux+Ux^2+Vx^2*3*Wx^2) * (Nx' * Nx) * w *le/2 ...
        + 2 * E * I *(Nxx' * Nxx)* w * le/2;
end

end