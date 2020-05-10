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
Nno=Nel+1;        % number of nodes in bend in linear_ansatz
% Nno=Nel*2+1;        % number of nodes in qudra-ansatz
le=l/Nel;         % length of an element
Iy=pi*(R^2-r^2)/4;  % Flaechentraegheitsmoment
Iz=Iy;

F=100;                      % force [N]
FVec= zeros(1*Nno,1);       % empty global force Vektor 
FVec(end)=F;


%% define empty matrice in u(x)
Kt=zeros(Nno);  % empty global stiffnes-matrix 
M=zeros(Nno);   % empty global mass-matrix 


Ae=zeros(2,1*Nno,Nel);
for ie=1:Nel
    Ae(1,1*(ie-1)+1,ie)=1;
    Ae(2,1*(ie-1)+2,ie)=1;
end

for j=1:5
    
   if j==1
       
       for k=1:Nel                                      % loop over every element
           
            %uka=uk0;
            %ue=Ae(:,:,j)*uka;                          % define diplacement-vector of current iteration
            ux=0;
            [Kte,Me] = Elementroutine_langs_linear(A,E,mu,rho,le,ux);
            Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
            M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
            %F =F+Ae(:,:,j)'*Fe
       end

       
   else
       
       for k=1:Nel                                     % loop over every element
           
           ux=(u(k+1)-u(k))/le;
           [Kte,Me] = Elementroutine_langs_linear(A,E,mu,rho,le,ux);
           Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
           M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
           
       end   
   end
   
   
   Kt(1,:) = [];
   Kt(:,1) = [];
   M(1,:) = [];
   M(:,1) = [];
   FVec(1) = [];

   u=Kt\FVec;
   u=[0;u];

   if j==5
       
   else
        Kt= zeros(1*Nno);                               % empty global stiffnes-matrix 
        M= zeros(1*Nno);                                % empty global mass-matrix 
        FVec= zeros(1*Nno,1);                           % empty global force Vektor 
        FVec(end)=F;
   end
   

     
end


%%
lVec = zeros (Nno, 1);     % vector with node coordinates

for k = 1: Nno
    lVec(k)=l/Nel*(k-1);
end

%% result by Solidworks
xls=xlsread('sta_mech_zug_test_mk01-only-Results-Displacement1-3');

%% plot

figure
grid on

subplot(1,2,1)
plot(lVec,u);
subplot(1,2,2)
plot(xls(:,2),xls(:,1));
    




