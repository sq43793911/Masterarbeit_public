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
Nel=1000;           % number of elements in bend
Nno=Nel+1;        % number of nodes in bend in linear_ansatz
% Nno=Nel*2+1;        % number of nodes in qudra-ansatz
le=l/Nel;         % length of an element
I=pi*(R^4-r^4)/4;  % Flaechentraegheitsmoment

Fx=100;                      % force [N]
Fy=100;                      % force [N]
Fz=100;                      % force [N]

FVec= zeros(5*Nel,1);       % empty global force Vektor 
FVec(end-4)=Fx;
FVec(end-3)=Fy;
FVec(end-1)=Fz;


uMat=[];
vMat=[];
vxMat=[];
wMat=[];
wxMat=[];


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
wx=zeros(Nno,1);
vx=zeros(Nno,1);

Nloop=5;

for j=1:Nloop
    
   if j==1
       
       for k=1:Nel                                      % loop over every element
            Ux=0;
            Vx=0;
            Wx=0;
            [Kte,Me] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I);
            Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
            M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
       end

       
   else
       
       for k=1:Nel                                     % loop over every element
           
           Ux=Fx/E/A;
           Vx=(vx(k+1)+vx(k))/2;
           Wx=(wx(k+1)+wx(k))/2;
           [Kte,Me] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I);
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
       v(m)=P(n+1);
       vx(m)=P(n+2);
       w(m)=P(n+3);
       wx(m)=P(n+4);
   end
       

   if j==Nloop
       
   else
        Kt= zeros(5*Nno);                               % empty global stiffnes-matrix 
        M= zeros(5*Nno);                                % empty global mass-matrix 
%         FVec= zeros(5*Nel,1);                           % empty global force Vektor 
%         FVec(end-4)=F;
%         FVec(end-3)=F;
%         FVec(end-1)=F;
   end
   
uMat=[uMat,u];
vMat=[vMat,v];
vxMat=[vxMat,vx];
wMat=[wMat,w];
wxMat=[wxMat,wx];
     
end


%%
lVec = zeros (Nno, 1);     % vector with node coordinates

for k = 1: Nno
    lVec(k)=l/Nel*(k-1);
end

%% result by Solidworks
u_xls_SW=xlsread('n.l.sta_mech_test_mk01-u.xlsx');
w_xls_SW=xlsread('n.l.sta_mech_test_mk01-w.xlsx');
v_xls_SW=xlsread('n.l.sta_mech_test_mk01-v.xlsx');

%% result by Ansys
u_xls_Ansys=xlsread('ansys_probe_u.xlsx');
w_xls_Ansys=xlsread('ansys_probe_w.xlsx');
v_xls_Ansys=xlsread('ansys_probe_v.xlsx');

%% plot

figure(1)
grid on

subplot(1,3,1)
plot(lVec,u);
subplot(1,3,2)
plot(lVec,w);
subplot(1,3,3)
plot(lVec,v);

figure(2)
subplot(1,3,1)
grid on
plot(u_xls_SW(:,3),u_xls_SW(:,2));
subplot(1,3,2)
plot(w_xls_SW(:,3),w_xls_SW(:,2));  
subplot(1,3,3)
plot(v_xls_SW(:,3),v_xls_SW(:,2));

figure(3)
grid on
subplot(1,3,1)
plot(u_xls_Ansys(:,3),u_xls_Ansys(:,5));
subplot(1,3,2)
plot(w_xls_Ansys(:,3),w_xls_Ansys(:,5));  
subplot(1,3,3)
plot(v_xls_Ansys(:,3),v_xls_Ansys(:,5));



