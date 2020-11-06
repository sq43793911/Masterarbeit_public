clear all
clc
close all

%###########################################
% main Program
%###########################################
% parameter
E=2.1e11;           % N/m^2
A=0.0001;           % m^2
l=10;               % m
rho=7850;           % Dichte in [kg/m^3]
mu=rho*A;           % Massenbelegung in [kg/m]
Nel=30;              % number of elements
Nno=Nel*3+1;          % number of nodes
le=l/Nel;           % length of an element



% define empty matrice
Kt=zeros(Nno);                                 % empty global stiffnes-matrix 
M=zeros(Nno);                                  % empty global mass-matrix 



% call element routinesmbclient
[Kte,Me] = Elementroutine_kubik_ohne_Abteilung(A,E,mu,le);

for j=1: 3 : Nno-2                                  % loop over every element
    
    M(j : j+3, j : j+3) = M(j : j+3, j : j+3) + Me;
    Kt(j : j+3, j : j+3) = Kt(j : j+3, j : j+3)+ Kte;
       
end

% implementation of essetial boundary conditions
Kt(1,:) = [  ];
Kt(:,1) = [  ];
M(1,:)  = [  ];
M(:,1)  = [  ];


% Kt(1,:) = [  ];
% Kt(:,1) = [  ];
% M(1,:)  = [  ];
% M(:,1)  = [  ];

% define system-matrix


null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins;
       -inv(M)*Kt,null];



% n = 2 * length (M);
% m = length (M);
% SysMat = zeros(n);
% x = size (SysMat, 1);
% SysMat(m+1: x ,1:Nel) = -inv(M)*Kt; 
% SysMat(1: m , m+1:x) = Eins;
% compute Eigenvalues
[V,D]=eig(SysMat);
temp_d = diag(D);
[nd, sortindex] = sort(temp_d);
temp_v = V(:,sortindex);
temp_f = imag(nd)/(2*pi);


% 
lamda = zeros(Nel, 1);
omega = zeros(Nel, 1);
f = zeros (Nel, 1);
for k = 1: Nel
    lamda (k) = (2 * k - 1 ) * pi /( 2 * l) ;
    omega (k) = lamda(k)*sqrt(E/rho);
    f(k) = omega(k)/(2*pi);
end