%###########################################
% Elementroutine
%###########################################
function [Kte,Me] = Elementroutine_kubik_ohne_Abteilung(A,E,mu,le)
% Elementroutine: compute Kte, Me
% define empty Kte
Kte=zeros(4);
 
% define empty Me 
Me=zeros(4);
 
%xiVec=[-sqrt(1/3),sqrt(1/3)];   % define sampling points for Gauss-quadrature      
%wVec =[1,1];   % weights for sampling points of Gauss-quadrature 


xiVec=[-sqrt(3/5),0,sqrt(3/5)];  % define sampling points for Gauss-quadrature
wVec=[5/9,8/9,5/9];             % weights for sampling points of Gauss-quadrature
for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
    % define N, B vector 
    N=[-1/16+xi/16+9*xi^2/16-9*xi^3/16   9/16-27*xi/16-9*xi^2/16+27*xi^3/16   9/16+27*xi/16-9*xi^2/16-27*xi^3/16   -1/16-xi/16+9*xi^2/16+9*xi^3/16];
    Nx=[1/16+9*xi/8-27*xi^2/16    -27/16-9*xi/8+81*xi^2/16   27/16-9*xi/8-81*xi^2/16   -1/16+9*xi/8+27*xi^2/16]*(2/le);

%     N=[-(1/4)+(5*xi)/8+(3*xi^2)/4-(9*xi^3)/8   3/5-(9*xi)/5-(3*xi^2)/5+(9*xi^3)/5   3/4+(9*xi)/8-(3*xi^2)/4-(9*xi^3)/8   -(1/10)+xi/20+(3*xi^2)/5+(9*xi^3)/20];
%     Nx=[5/8+(3*xi)/2-(27*xi^2)/8    (27*xi^2)/5-(6*xi)/5-9/5   9/8-(3*xi)/2-(27*xi^2)/8   1/20+(6*xi)/5+(27*xi^2)/20]*(2/le);

% 
% %      
%     N = [-0.16666666666666666 + 0.16666666666666666*xi + ... 
%     0.6666666666666666*xi^2 - 0.6666666666666667* ...
%       xi^3    0.6666666666666666 - 1.3333333333333333*xi - ...
%     0.6666666666666666*xi^2 + 1.3333333333333333*...
%       xi^3   0.6666666666666666 + 1.3333333333333333*xi - ...
%     0.6666666666666666*xi^2 - 1.3333333333333333*...
%       xi^3    -0.16666666666666666 - 0.16666666666666669*xi + ...
%     0.6666666666666666*xi^2 + 0.6666666666666667*...
%       xi^3];
%     Nx = [0.16666666666666666 + 1.3333333333333333*xi - ...
%      2.*xi^2   -1.3333333333333333 - 1.3333333333333333*xi + ...
%      4.*xi^2   1.3333333333333333 - 1.3333333333333333*xi - ...
%      4.*xi^2   -0.16666666666666669 + 1.3333333333333333*xi + ...
%      2.*xi^2]*(2/le);

    % compute Kte and Me for sampling point of Gauss-integration
    Me=Me + mu * (N' * N) * w;
    Kte=Kte + E * A * (Nx' * Nx) * w;
end

end