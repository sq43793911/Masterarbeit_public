%###########################################
% Elementroutine
%###########################################
function [Ktue,Ktwe,Ktve,Mue,Mwe,Mve] = Elementroutine_Hohr(A,E,mu,rho,le,Iy,Iz)
% Elementroutine: compute Kte, Me
% define empty Kte
Ktue=zeros(3);
Ktwe=zeros(4);
Ktve=zeros(4);
Kte=zeros(4);
% define empty Me 
Mue=zeros(3);
Mwe=zeros(4);
Mve=zeros(4);
Me=zeros(4);
%xiVec=[-sqrt(1/3),sqrt(1/3)];   % define sampling points for Gauss-quadrature      
%wVec =[1,1];   % weights for sampling points of Gauss-quadrature 


% xiVec=[-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];  % define sampling points for Gauss-quadrature
% wVec=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];   % weights for sampling points of Gauss-quadrature

xiVec=[-(1/3)*sqrt(5+2*(sqrt(10/7))), -(1/3)*sqrt(5-2*sqrt(10/7)), 0, (1/3)*sqrt(5-2*sqrt(10/7)), (1/3)*sqrt(5+2*sqrt(10/7)) ];
wVec=[(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900 ];

for i=1:length(xiVec)
    xi=xiVec(i);
    w =wVec(i);
    
    %% define N, B vector von langsverschiebung
%     Nu=[-1/16+xi/16+9*xi^2/16-9*xi^3/16   9/16-27*xi/16-9*xi^2/16+27*xi^3/16   9/16+27*xi/16-9*xi^2/16-27*xi^3/16   -1/16-xi/16+9*xi^2/16+9*xi^3/16];
%     Nux=[1/16+9*xi/8-27*xi^2/16    -27/16-9*xi/8+81*xi^2/16   27/16-9*xi/8-81*xi^2/16   -1/16+9*xi/8+27*xi^2/16]*(2/le);

    Nu=[xi^2/2-xi/2  1-xi^2  xi/2+xi^2/2];
    Nux=[xi-0.5  -2*xi  0.5+xi]*(2/le);
    % compute Kute and Mue for sampling point of Gauss-integration
    Mue=Mue + mu * (Nu' * Nu) * w;
    Ktue=Ktue + E * A * (Nux' * Nux) * w;
    
    %% define N, B vector von y
    Nw = [1/2-(3*xi)/4+(xi^3)/4   1/4-xi/4-(xi^2)/4+(xi^3)/4   1/2+(3*xi)/4-(xi^3)/4   -1/4-xi/4+(xi^2)/4+(xi^3)/4];
    Nwxx = [(3*xi)/2   -1/2+(3*xi)/2   -(3*xi)/2   1/2+(3*xi)/2]*((2/le)^2);
     
    % compute Kwte and Mwe for sampling point of Gauss-integration
    Mwe=Me + rho * A * (Nw' * Nw) * w;
    Ktwe=Kte + E * Iy * (Nwxx' * Nwxx) * w;
    
    %% define N, B vector von z
    Nv = [1/2-(3*xi)/4+(xi^3)/4   1/4-xi/4-(xi^2)/4+(xi^3)/4   1/2+(3*xi)/4-(xi^3)/4   -1/4-xi/4+(xi^2)/4+(xi^3)/4];
    Nvxx = [(3*xi)/2   -1/2+(3*xi)/2   -(3*xi)/2   1/2+(3*xi)/2]*((2/le)^2);
     
    % compute Kvte and Mve for sampling point of Gauss-integration
    Mve=Me + rho * A * (Nv' * Nv) * w;
    Ktve=Kte + E * Iz * (Nvxx' * Nvxx) * w;
    
%     %% compute Kte and Me 
%     Me=Mue+Mwe+Mve;
%     Kte=Ktue+Ktwe+Ktve;
    
end

end