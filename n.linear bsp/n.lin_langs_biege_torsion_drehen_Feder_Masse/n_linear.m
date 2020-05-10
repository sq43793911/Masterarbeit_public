clear
close all
%###########################################
% main Program
%###########################################
%% parameter
E=2.1e11;         % N/m^2
D0=0.01;           % Durchmesser m
R=D0/2;
d=0.002;          % Wandstaeker m
r=R-d;
A=pi*(R^2-r^2);   % Flaeche m^2
l=0.34;           % m
rho=7800;         % Dichte in [kg/m^3]
mu=rho*A;         % Massenbelegung in [kg/m]
Nel=100;           % number of elements in bend
Nno=Nel+1;        % number of nodes in bend in linear_ansatz
% Nno=Nel*2+1;        % number of nodes in qudra-ansatz
le=l/Nel;         % length of an element
I=pi*(R^4-r^4)/4;  % Flaechentraegheitsmoment
It=2*I;     % Torsionstraegheitsmoment
v=0.28;                % Poissonzahl
G=E/(2*(1+v));        % Schubmodul
c=0;                  % Feder-Steifigkeit [N/m]
m=0.0145;                % Masse  [kg]

Omega=60;    % Drehgeschwindigkeit [r/s]
e=-0.002;      % Exzentrizitaet      [m]
q=6;          % Freiheitsgrad

Fx=0;                      % force [N]
Fy=0;                      % force [N]
Fz=0;                      % force [N]
M=0;                       % moment [N*m]

FVec= zeros(q*Nel,1);       % empty global force Vektor 
FVec(end-5)=Fx;
FVec(end-4)=Fy;
FVec(end-2)=Fz;
FVec(end)=M;


uMat=[];
vMat=[];
vxMat=[];
wMat=[];
wxMat=[];


%% define empty matrice in u(x)
Kt=zeros(Nno*q);  % empty global stiffnes-matrix 
M=zeros(Nno*q);   % empty global mass-matrix 
B=zeros(Nno*q);   
Q=zeros(Nno*q);


Ae=zeros(2*q,q*Nno,Nel);
for ie=1:Nel
    for i=1:2*q
    Ae(i,q*(ie-1)+i,ie)=1;
    end
end

u=zeros(Nno,1);
w=zeros(Nno,1);
v=zeros(Nno,1);
wx=zeros(Nno,1);
vx=zeros(Nno,1);
phi=zeros(Nno,1);

Nloop=10;
Ux=0;
Vx=0;
Wx=0;

for j=1:Nloop
    
   if j==1
       
       for k=1:Nel                                      % loop over every element
            [Kte,Me,Be,Qe,Ce,MeM,KteM,BeM,QeM] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I,It,G,Omega,e,c,m);
            if k==Nel
                Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k)+Ae(:,:,k)'*Ce*Ae(:,:,k)+Ae(:,:,k)'*KteM*Ae(:,:,k);
                M =M+Ae(:,:,k)'*Me*Ae(:,:,k)+Ae(:,:,k)'*MeM*Ae(:,:,k);              
                B =B+Ae(:,:,k)'*Be*Ae(:,:,k)+Ae(:,:,k)'*BeM*Ae(:,:,k);
                Q =Q+Ae(:,:,k)'*Qe*Ae(:,:,k)+Ae(:,:,k)'*QeM*Ae(:,:,k);
            else
            Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
            M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
            B =B+Ae(:,:,k)'*Be*Ae(:,:,k);
            Q =Q+Ae(:,:,k)'*Qe*Ae(:,:,k);
            end
       end

       
   else
       
       for k=1:Nel                                     % loop over every element
           
%                       Ux=Fx/E/A;
%                       Vx=(vx(k+1)+vx(k))/2;
%                       Wx=(wx(k+1)+wx(k))/2;
           
           Ux=(u(k+1)-u(k))/le;
           Wx=(w(k+1)-w(k))/le;
           Vx=(v(k+1)-v(k))/le;
           [Kte,Me,Be,Qe,Ce,MeM,KteM,BeM,QeM] = Elementroutine_n_linear(A,E,rho,le,Ux,Vx,Wx,I,It,G,Omega,e,c,m);
           if k==Nel
               Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k)+Ae(:,:,k)'*Ce*Ae(:,:,k)+Ae(:,:,k)'*KteM*Ae(:,:,k);
               M =M+Ae(:,:,k)'*Me*Ae(:,:,k)+Ae(:,:,k)'*MeM*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
               B =B+Ae(:,:,k)'*Be*Ae(:,:,k)+Ae(:,:,k)'*BeM*Ae(:,:,k);
               Q =Q+Ae(:,:,k)'*Qe*Ae(:,:,k)+Ae(:,:,k)'*QeM*Ae(:,:,k);
           else
               Kt=Kt+Ae(:,:,k)'*Kte*Ae(:,:,k);             % place the distribution of every element to right place in global stiffnes-matrix
               M =M+Ae(:,:,k)'*Me*Ae(:,:,k);               % place the distribution of every element to right place in global mass-matrix
               B =B+Ae(:,:,k)'*Be*Ae(:,:,k);
               Q =Q+Ae(:,:,k)'*Qe*Ae(:,:,k);
           end
       end
   end
   
   for zz=1:q
       Kt(1,:) = [];
       Kt(:,1) = [];
       M(1,:) = [];
       M(:,1) = [];
       B(1,:) = [];
       B(:,1) = [];
   end
   
   
   QVec=diag(Q);
   QVec(1:q)=[];
   
   AVec=FVec+QVec;

   P=Kt\AVec;
   for yy=1:q
       P=[0;P];
   end
   
   
   
   for xx=1:Nno
       n=(xx-1)*q+1;
       u(xx)=P(n);
       w(xx)=P(n+1);
       wx(xx)=P(n+2);
       v(xx)=P(n+3);
       vx(xx)=P(n+4);
       phi(xx)=P(n+5);
   end
       

   if j==Nloop
       
   else
        Kt= zeros(q*Nno);                               % empty global stiffnes-matrix 
        M = zeros(q*Nno);                                % empty global mass-matrix 
        B = zeros(q*Nno);   
        Q = zeros(q*Nno);
   end
   
% uMat=[uMat,u];
% vMat=[vMat,v];
% vxMat=[vxMat,vx];
% wMat=[wMat,w];
% wxMat=[wxMat,wx];
     
end

%%
% define system-matrix
null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins; -inv(M)*Kt,-inv(M)*B];

[V,D]=eig(SysMat);

temp_d = diag(D);
org_d = temp_d;
for i=1:Nel*6
 temp_d(i,:)=[];
end

[nd, sortindex] = sort(temp_d);
temp_v = V(:,sortindex);
temp_f = -imag(nd)/(2*pi);

%%
lVec = zeros (Nno, 1);     % vector with node coordinates

for k = 1: Nno
    lVec(k)=l/Nel*(k-1);
end


%% plot

figure('Name','Programm')
grid on

subplot(1,3,1)
plot(lVec,u);
title('u(x)')
subplot(1,3,2)
plot(lVec,w);
title('w(x)')
subplot(1,3,3)
plot(lVec,v);
title('v(x)')

%% result by Solidworks
% u_xls_SW=xlsread('n.l.sta_mech_test_mk01-u.xlsx');
% w_xls_SW=xlsread('n.l.sta_mech_test_mk01-w.xlsx');
% v_xls_SW=xlsread('n.l.sta_mech_test_mk01-v.xlsx');
% 
% figure('Name','Solidworks')
% subplot(1,3,1)
% grid on
% plot(u_xls_SW(:,3),u_xls_SW(:,2));
% title('u(x)')
% subplot(1,3,2)
% plot(w_xls_SW(:,3),w_xls_SW(:,2));  
% title('w(x)')
% subplot(1,3,3)
% plot(v_xls_SW(:,3),v_xls_SW(:,2));
% title('v(x)')


%% result by Ansys
% u_xls_Ansys=xlsread('ansys_probe_u.xlsx');
% w_xls_Ansys=xlsread('ansys_probe_w.xlsx');
% v_xls_Ansys=xlsread('ansys_probe_v.xlsx');
% figure('Name','Ansys')
% grid on
% subplot(1,3,1)
% plot(u_xls_Ansys(:,3),u_xls_Ansys(:,5));
% title('u(x)')
% subplot(1,3,2)
% plot(w_xls_Ansys(:,3),w_xls_Ansys(:,5));
% title('w(x)')
% subplot(1,3,3)
% plot(v_xls_Ansys(:,3),v_xls_Ansys(:,5));
% title('v(x)')


% u_xls_Ansys=xlsread('ansys_probe_u_new.xlsx');
% w_xls_Ansys=xlsread('ansys_probe_w_new.xlsx');
% v_xls_Ansys=xlsread('ansys_probe_v_new.xlsx');
% figure('Name','Ansys')
% grid on
% subplot(1,3,1)
% plot(u_xls_Ansys(:,3),u_xls_Ansys(:,5));
% title('u(x)')
% subplot(1,3,2)
% plot(w_xls_Ansys(:,3),w_xls_Ansys(:,5));
% title('w(x)')
% subplot(1,3,3)
% plot(v_xls_Ansys(:,3),v_xls_Ansys(:,5));
% title('v(x)')

