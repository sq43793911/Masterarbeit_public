%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_quadra(A,E,mu,le)
% Elementroutine: compute Kte, Me
% define empty Kte
Kte=zeros(3); 
% define empty Me 
Me=zeros(3);

% define sampling points for Gauss-quadrature
xiVec=[-sqrt(1/3),sqrt(1/3)]; 
% weights for sampling points of Gauss-quadrature 
wVec =[1,1];   

for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
    % define N, B vector 
    N=[xi^2/2-xi/2  1-xi^2  xi/2+xi^2/2];
    Nx=[xi-0.5  -2*xi  0.5+xi]*(2/le);
        
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + mu * (N' * N) * w;
    Kte=Kte + E * A * (Nx' * Nx) * w;
end

end