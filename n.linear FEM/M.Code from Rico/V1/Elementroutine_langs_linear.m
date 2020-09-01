%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_langs_linear(A,E,mu,rho,le,ux)
% Elementroutine: compute Kte, Me
%%
%define empty Kte
Kte=[0,0;
     0,0]; 
% define empty Me 
Me=[0,0;
    0,0];
%%
% % define empty Kte qudra
% Kte=zeros(3); 
% % define empty Me qudra
% Me=zeros(3);

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
    N=[0.5-xi/2  0.5+xi/2 ];
    Nx=[-0.5  0.5]*(2/le);
    
    % define N, B vector qudra
%     N=[xi^2/2-xi/2  1-xi^2  xi/2+xi^2/2];
%     Nx=[xi-0.5  -2*xi  0.5+xi]*(2/le);
        
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + rho * A * (N' * N) * w*le/2;
    Kte=Kte + E * A * (1+(3*ux^2)/2+3*ux) * (Nx' * Nx) * w*le/2;
end

end