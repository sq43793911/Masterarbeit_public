%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_linear(A,E,rho,le)
% Elementroutine: compute Kte, Me
% define empty Kte
Kte=[0,0;
     0,0]; 
% define empty Me 
Me=[0,0;
    0,0];
 
xiVec=[-sqrt(1/3),sqrt(1/3)];   % define sampling points for Gauss-quadrature      
wVec =[1,1];   % weights for sampling points of Gauss-quadrature 

for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
    % define N, B vector 
    N=[0.5-xi/2  0.5+xi/2 ];
    Nx=[-0.5  0.5]*(2/le);
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + rho * A * N' * N * w;
    Kte=Kte + E * A * Nx' * Nx * w;
end

end