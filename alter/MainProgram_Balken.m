clear
clc
close all
%###########################################
% main Program
%###########################################
% parameter
E=2.1e11;        % N/m^2
A=0.0001;        % m^2
l=10;            % m
rho=7850;        % Dichte in [kg/m^3]
mu=rho*A;        % Massenbelegung in [kg/m]
Nel=30;         % number of elements
Nno=Nel+1;       % number of nodes
le=l/Nel;        % length of an element
B=0.005;         % Breit von Balken
H=0.02;          % Hoehe von Balken
I=(B*(H^3))/12;  % Flaechentraegheitsmoment

% define empty matrice
Kt=zeros(Nno*2);  % empty global stiffnes-matrix 
M=zeros(Nno*2);   % empty global mass-matrix 

% call element routinesmbclient
[Kte,Me] = Elementroutine_Balken(A,E,rho,le,I);

% loop over every element
for j=1: 2 : length(M)-2                                   
    M(j : j+3, j : j+3) = M(j : j+3, j : j+3) + Me;
    Kt(j : j+3, j : j+3) = Kt(j : j+3, j : j+3)+ Kte;
end

% implementation of essetial boundary conditions
Kt(1,:) = [  ];
Kt(:,1) = [  ];
M(1,:)  = [  ];
M(:,1)  = [  ];

Kt(1,:) = [  ];
Kt(:,1) = [  ];
M(1,:)  = [  ];
M(:,1)  = [  ];

% define system-matrix
null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins; -inv(M)*Kt,null];

% compute Eigenvalues
[V,D]=eig(SysMat);

temp_d = diag(D);
[nd, sortindex] = sort(temp_d);
temp_v = V(:,sortindex);
temp_f = imag(nd)/(2*pi);

for i=1:Nel*2
 temp_f(i,:)=[];
end

lVec = zeros (Nno, 1);     % vector with node coordinates

for k = 1: Nno
    lVec(k)=l/Nel*(k-1);
end

EigMat=temp_v(1:length(M),:) ;   % delete lover part von eigenvectors 
EigMat=imag(EigMat);
for k = 1: length(M)/2      % delete unnecessary rows 
   EigMat(1+k,:)=[];       
end

for k = 1: length(M)
   EigMat(:,k)=[];         % delete unnecessary colums
end

EigMat=[zeros(1,length(M));EigMat];  % fill up EigMat for first with zero displacement

figure
hold on
grid on
plot(lVec,-EigMat(:,1),'LineWidth',4);
plot(lVec,EigMat(:,2),'--','LineWidth',4);
plot(lVec,EigMat(:,3),':','LineWidth',4);
set(gca,'FontSize',24);
legend('1-Mode', '2-Mode', '3-Mode','Location','northwest');

% analysis method
lamda = zeros(Nel, 1);
omega = zeros(Nel, 1);
f = zeros (Nel, 1);
for k = 1: Nel
    switch k
        case 1 
            lamda(k) = 1.87510 ;
        case 2
            lamda(k) = 4.69409 ;
        case 3
            lamda(k) = 7.85476 ;
        case 4
            lamda(k) = 10.99554;
    end
    if k >= 5
            lamda(k) = ((2*k-1)*pi)/2 ;
    end
      
    omega (k) = ((lamda(k))^2)*sqrt((E*I)/(rho*A*(l^4)));
    f(k)=omega(k)/(2*pi);
end

