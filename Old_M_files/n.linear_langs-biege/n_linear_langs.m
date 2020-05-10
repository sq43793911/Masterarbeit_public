clear
clc
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
Nno=Nel+1;        % number of nodes in bend in linear_ansatz
% Nno=Nel*2+1;        % number of nodes in qudra-ansatz
le=l/Nel;         % length of an element
I=pi*(R^2-r^2)/4;  % Flaechentraegheitsmoment

F=100;                      % force [N]
FVec= zeros(5*Nel,1);       % empty global force Vektor 
FVec(end-4)=F;
FVec(end-3)=F;
FVec(end-1)=F;


%% define empty matrice in u(x)
Kt=zeros(Nno*5);  % empty global stiffnes-matrix 
M=zeros(Nno*5);   % empty global mass-matrix 


Ae=zeros(10,5*Nno,Nel);
for ie=1:Nel
    for i=1:10
    Ae(i,5*(ie-1)+i,ie)=1;
    end
end

u=zeros(Nno,1);
w=zeros(Nno,1);
v=zeros(Nno,1);

for j=1:5
    
   if j==1
       
       for k=1:Nel                                      % loop over every element
            Ux=0;
            Vx=0;
            Wx=0;
            [Kte,Me] = Elementroutine_langs_n_linear(A,E,rho,le,Ux,Vx,Wx,I);
            Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
            M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
       end

       
   else
       
       for k=1:Nel                                     % loop over every element
           
           Ux=(u(k+1)-u(k))/le;
           Wx=(w(k+1)-w(k))/le;
           Vx=(v(k+1)-v(k))/le;
           [Kte,Me] = Elementroutine_langs_n_linear(A,E,rho,le,Ux,Vx,Wx,I);
           Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
           M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
           
       end   
   end
   
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];

   P=Kt\FVec;
   P=[0;0;0;0;0;P];
   
   
   for m=1:Nno
       n=(m-1)*5+1;
       u(m)=P(n);
       w(m)=P(n+1);
       v(m)=P(n+3);
   end
       

   if j==5
       
   else
        Kt= zeros(5*Nno);                               % empty global stiffnes-matrix 
        M= zeros(5*Nno);                                % empty global mass-matrix 
        FVec= zeros(5*Nel,1);                           % empty global force Vektor 
        FVec(end-4)=F;
        FVec(end-3)=F;
        FVec(end-1)=F;
   end
   

     
end


%%
lVec = zeros (Nno, 1);     % vector with node coordinates

for k = 1: Nno
    lVec(k)=l/Nel*(k-1);
end

%% result by Solidworks
% xls=xlsread('sta_mech_zug_biege_test_mk01-w');

%% plot

figure
grid on

subplot(1,3,1)
plot(lVec,u);
subplot(1,3,2)
plot(lVec,w);
subplot(1,3,3)
plot(lVec,v);
% subplot(1,4,4)
% plot(xls(:,3),xls(:,2));
    




