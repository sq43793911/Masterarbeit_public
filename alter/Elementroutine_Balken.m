%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_Balken(A,E,rho,le,I)
% Elementroutine: compute Kte, Me
% define empty Kte
Kte=zeros(4);
% define empty Me 
Me=zeros(4);
% define sampling points for Gauss-quadrature
xiVec=[-sqrt(3/5),0,sqrt(3/5)]; 
% weights for sampling points of Gauss-quadrature
wVec=[5/9,8/9,5/9];            
for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
    % zwei Freiheitsgrade
    N = [1/2-(3*xi)/4+(xi^3)/4   1/4-xi/4-(xi^2)/4+(xi^3)/4   1/2+(3*xi)/4-(xi^3)/4   -1/4-xi/4+(xi^2)/4+(xi^3)/4];
    Nxx = [(3*xi)/2   -1/2+(3*xi)/2   -(3*xi)/2   1/2+(3*xi)/2]*((2/le)^2);
     
    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + rho * A * (N' * N) * w;
    Kte=Kte + E * I * (Nxx' * Nxx) * w;
end

end