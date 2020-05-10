clear
close
clc
%###########################################
% main Program
%###########################################
% parameter
E=2.1e11;           % N/m^2
A=0.0001;           % m^2
l=10;               % m
rho=7850;           % Dichte in [kg/m^3]
mu=rho*A;           % Massenbelegung in [kg/m]
Nel=30;             % number of elements
Nno=Nel*2+1;        % number of nodes
le=l/Nel;           % length of an element

% define empty matrice
Kt=zeros(Nno);   % empty global stiffnes-matrix 
M=zeros(Nno);    % empty global mass-matrix 



% call element routinesmbclient
[Kte,Me] = Elementroutine_quadra(A,E,mu,le);
% loop over every element
for j=1:2 : Nno-1                                    
    M(j : j+2, j : j+2) = M(j : j+2, j : j+2) + Me;
    Kt(j : j+2, j : j+2) = Kt(j : j+2, j : j+2)+ Kte;   
end

% implementation of essetial boundary conditions
Kt(1,:) = [  ];
Kt(:,1) = [  ];
M(1,:)  = [  ];
M(:,1)  = [  ];

% define system-matrix
null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins;
        -inv(M)*Kt,null];
   
% compute Eigenvalues
[V,D]=eig(SysMat);
temp_d = diag(D);
[nd, sortindex] = sort(temp_d);
temp_v = V(:,sortindex);
temp_f = imag(nd)/(2*pi);

% analysis method
lamda = zeros(18, 1);
omega = zeros(18, 1);
f = zeros (18, 1);
for k = 1: 18
    lamda (k) = (2 * k - 1 ) * pi /( 2 * l) ;
    omega (k) = lamda(k)*sqrt(E/rho);
    f(k) = omega(k)/(2*pi);
end