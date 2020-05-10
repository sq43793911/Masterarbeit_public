clear
clc
close
%###########################################
% main Program
%###########################################
%% parameter
E=2.1e11;         % N/m^2
D=0.01;           % Durchmesser m
R=D/2;
d=0.002;          % Wandstaeker m
r=(D-d)/2;
A=pi*(R^2-r^2);   % Flaeche m^2
l=0.27;           % m
rho=7850;         % Dichte in [kg/m^3]
mu=rho*A;         % Massenbelegung in [kg/m]
Nel=20;           % number of elements in bend
Nno=Nel+1;        % number of nodes in bend
le=l/Nel;         % length of an element
Iy=pi*(R^2-r^2)/4;  % Flaechentraegheitsmoment
Iz=Iy;

Nuel=30;             % number of elements in u(x) 
Nuno=Nel*2+1;        % number of nodes in u(x)

%% define empty matrice in u(x)
Ktu=zeros(Nuno);  % empty global stiffnes-matrix 
Mu=zeros(Nuno);   % empty global mass-matrix 

% call element routinesmbclient
[Ktue,Ktwe,Ktve,Mue,Mwe,Mve] = Elementroutine_Hohr(A,E,mu,rho,le,Iy,Iz);

% loop over every element
for j=1:2 : Nno-1                                    
    Mu(j : j+2, j : j+2) = Mu(j : j+2, j : j+2) + Mue;
    Ktu(j : j+2, j : j+2) = Ktu(j : j+2, j : j+2)+ Ktue;   
end

%% define empty matrice in bend
Ktw=zeros(Nno*2);  % empty global stiffnes-matrix 
Mw=zeros(Nno*2);   % empty global mass-matrix 

% loop over every element
for j=1: 2 : length(Mw)-2                                   
    Mw(j : j+3, j : j+3) = Mw(j : j+3, j : j+3) + Mwe;
    Ktw(j : j+3, j : j+3) = Ktw(j : j+3, j : j+3)+ Ktwe;
end

Ktv=zeros(Nno*2);  % empty global stiffnes-matrix 
Mv=zeros(Nno*2);   % empty global mass-matrix 

% loop over every element
for j=1: 2 : length(Mv)-2                                   
    Mv(j : j+3, j : j+3) = Mv(j : j+3, j : j+3) + Mve;
    Ktv(j : j+3, j : j+3) = Ktv(j : j+3, j : j+3)+ Ktve;
end

%% implementation of essetial boundary conditions

Ktu(1,:) = [  ];
Ktu(:,1) = [  ];
Mu(1,:)  = [  ];
Mu(:,1)  = [  ];

Ktw(1,:) = [  ];
Ktw(:,1) = [  ];
Mw(1,:)  = [  ];
Mw(:,1)  = [  ];

Ktw(1,:) = [  ];
Ktw(:,1) = [  ];
Mw(1,:)  = [  ];
Mw(:,1)  = [  ];

Ktv(1,:) = [  ];
Ktv(:,1) = [  ];
Mv(1,:)  = [  ];
Mv(:,1)  = [  ];

Ktv(1,:) = [  ];
Ktv(:,1) = [  ];
Mv(1,:)  = [  ];
Mv(:,1)  = [  ];

%% amount all matrix

M=Mu+Mw+Mv;
Kt=Ktu+Ktw+Ktv;

%% define system-matrix
null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins; 
       -inv(M)*Kt,null];

%% compute Eigenvalues
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

%% analysis method
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
      
    omega (k) = ((lamda(k))^2)*sqrt((E*Iy)/(rho*A*(l^4)));
    f(k)=omega(k)/(2*pi);
end

