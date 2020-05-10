clear
clc
close
%###########################################
% main Program
%###########################################
% parameter
E=2.1e11;       % N/m^2
A=0.0001;       % m^2
L=2;            % m  length of plate
B=2;            % m  wide of plate
rho=7850;       % Dichte in [kg/m^3]
mu=rho*A;       % Massenbelegung in [kg/m]
Nelx=15;        % number of elements per rand
Nely=15;
Nnox=Nelx+1;    % number of nodes per rand
Nnoy=Nely+1;
Nel_all=Nelx * Nely;  % total elements
Nno_all=Nnox * Nnoy;  % total nodes
lex=L/Nelx;           % length of an element
ley=B/Nely;
h=0.01;          % m Dicke von Platte
v=0.30;               % Poissionzahl

% define empty matrice
Kt=zeros(Nno_all*4);     % empty global stiffnes-matrix 
M=zeros(Nno_all*4);      % empty global mass-matrix 

% call element routinesmbclient
[Kte,Me] = Elementroutine_Platten(h,E,rho,lex,ley,v);
% loop over every element
r=0;
a=zeros(4,1);
b=a;
c=[1 5 9 13];
d=[4 8 12 16];
for j=1 : Nely                                  
    for k=1:Nelx
        a(1)=((j-1)*Nnox+k-1)*4+1;
        b(1)=a(1)+3;
        a(2)=((j-1)*Nnox+k)*4+1;
        b(2)=a(2)+3;
        a(3)=((j*Nnox)+k-1)*4+1;
        b(3)=a(3)+3;
        a(4)=((j*Nnox)+k)*4+1;
        b(4)=a(4)+3;
        for m=1:4
            for n=1:4
                M(a(m):b(m), a(n):b(n)) = M(a(m):b(m), a(n):b(n)) + Me(c(m):d(m),c(n):d(n));
                Kt(a(m):b(m), a(n):b(n)) = Kt(a(m):b(m), a(n):b(n)) + Kte(c(m):d(m),c(n):d(n));
            end
        end
    end
end

% implementation of essetial boundary conditions
M(1:Nnox*4,:) = [];
M(:,1:Nnox*4) = [];

Kt(1:Nnox*4,:) = [];
Kt(:,1:Nnox*4) = [];

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

%% Mode Plot
% vector with node coordinates
xVec = zeros (Nnox, 1);     
yVec = zeros (Nnoy, 1); 
for k = 1: Nnox
    xVec(k)=L/Nelx*(k-1);    
end

for g = 1: Nnoy
    yVec(g)=B/Nely*(g-1);    
end

% delete lover part von eigenvectors 
EigMat=temp_v(1:length(M),:) ;   
EigMat=imag(EigMat);
% delete unnecessary rows 
for k = 1: length(M)/4      
   EigMat(1+k:1+k+2,:)=[];       
end
% delete unnecessary colums
for k = 1: length(M)
   EigMat(:,k)=[];         
end
% fill up EigMat
p = zeros(Nnox,Nnoy-1,9);
for i = 1:9
    for j = 1:Nnoy-1
        for k = 1:Nnox
            p(k,j,i) = EigMat(k+(j-1)*Nnox,i);
        end
    end
end

for i = 1 : 9
    p1(:,:,i) = [zeros(Nnox,1),p(:,:,i)];
end

figure
hold on
for k = 1:9
    subplot(3,3,k);
    surf(yVec,xVec,p1(:,:,k));
    set(gca,'FontSize',20);
    title([num2str(k),'-Mode']);
end