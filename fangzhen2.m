clc
clear

N=100;
slid1=zeros(1,N+1);
s1bar=zeros(1,N+1);

A=[0.7248    0.1752
        -1.9131    0.3744];
B=[0.0252
        0.1752];
g1 = [-0.2520 -0.0323];
g2 = [-0.6345,0.1587];
g3 = [-0.1587,0.1603];
g4 = [-0.1587,0.1603];
phi=[2,0;0,4];
phi0=[0.2,0;0,0.3];
K1=[0.0108 -0.1047];
K2=[1.0194 -0.0534];
K3=[1.0194 -0.0534];
K4=[1.0194 -0.0534];
miu=0.03;
Emax=1.5;
Emin=0.5;

X0ini=[2;2];
X0state=zeros(2,N);
X0=X0state;              
X0(:,1)=X0ini;
X0set=zeros(2,N);

X1ini=[8;-6];
X1state=zeros(2,N);
X1=X1state;              
X1(:,1)=X1ini;
x1=X1ini;
last_X1=X1ini;
Q1=zeros(1,N);

T0_1=zeros(1,N);
T0_2=zeros(1,N);

T1_1=zeros(1,N);
T1_2=zeros(1,N);

T2_1=zeros(1,N);
T2_2=zeros(1,N);

T3_1=zeros(1,N);
T3_2=zeros(1,N);


R1_1=zeros(1,N);
R1_2=zeros(1,N);

R2_1=zeros(1,N);
R2_2=zeros(1,N);

R3_1=zeros(1,N);
R3_2=zeros(1,N);

X2ini=[-2;1];
X2state=zeros(2,N);
X2=X2state;              
X2(:,1)=X2ini;
x2=X2ini;
last_X2=X2ini;
Q2=zeros(1,N);

X3ini=[-5;3];
X3state=zeros(2,N);
X3=X3state;              
X3(:,1)=X3ini;
x3=X3ini;
last_X3=X3ini;
Q3=zeros(1,N);

X4ini=[1;3];
X4state=zeros(2,N);
X4=X4state;              
X4(:,1)=X4ini;
x4=X4ini;
last_X4=X4ini;
Q4=zeros(1,N);

Eta1ini=0;
Eta1state=zeros(1,N);
Eta1=Eta1state;
Eta1(:,1)=Eta1ini;
eta1=Eta1ini;
last_Eta1=Eta1ini;

Eta2ini=0;
Etastate2=zeros(1,N);
Eta2=Etastate2;
Eta2(:,1)=Eta2ini;
eta2=Eta2ini;
last_Eta2=Eta2ini;

Eta3ini=0;
Etastate3=zeros(1,N);
Eta3=Etastate3;
Eta3(:,1)=Eta3ini;
eta3=Eta3ini;
last_Eta3=Eta3ini;

u0ini=0;
u0state=zeros(1,N);
u0=u0state;              
u0(:,1)=u0ini;
u1ini=0;
u1state=zeros(1,N);
u1set=u1state;              
u1(:,1)=u1ini;

fini=0;
fstate=zeros(1,N);
f=fstate;              
f(:,1)=fini;

juzhen1=zeros(2,2);
juzhen2=zeros(2,2);
juzhen3=zeros(2,2);



a=[1,0;0,0];
aa=[0,0;0,1];

    for k=1:N 
  
   u0(:,k)=0.2*cos(k*pi);
   f(:,k)=-0.04*sin(((k-1)/(k+1)));
   X0(:,k+1)=A*X0(:,k)+B*(u0(:,k)+f(:,k));
   

       if mod(k,2)==1
       barcer1=aa*(X1(:,k)-X0(:,k));
       celiang1=aa*(X1(:,k)-last_X1);
       else
       barcer1=a*(X1(:,k)-X0(:,k));
       celiang1=a*(X1(:,k)-last_X1);
       end
       
       s1bar(:,k)=g1*(barcer1);
       u1(:,k)=-K1*(barcer1)+u0(:,k)-miu*norm(barcer1)*s1bar(:,k)/(norm(s1bar(:,k))+0.0005)
       
       X1(:,k+1)=A*X1(:,k)+B*u1(:,k);      
       Eta1(:,1)=0.5;
       Eta1(:,k+1)=0.01*Eta1(:,k)+(Emax*barcer1'*phi*barcer1-(celiang1)'*phi0*(celiang1));
       
       if (celiang1)'*phi0*(celiang1)>Emax*barcer1'*phi*barcer1+0*Eta1(:,k)
           Q1(:,k)=1;
           T1_1(:,k)=1;
           T0_1(:,k)=-2;
           if (aa*(X1(:,k)-last_X1))'*(aa*(X1(:,k)-last_X1))>=(a*(X1(:,k)-last_X1))'*(a*(X1(:,k)-last_X1))
               R1_2(:,k)=2;
               juzhen1=aa;
           else
               R1_1(:,k)=1;
               juzhen1=a;  
           end
           last_X1=X1(:,k);

       elseif(celiang1)'*phi0*(celiang1)>=Emin*barcer1'*phi*barcer1
           last_X1=X1(:,k);
           Q1(:,k)=1;
           T1_2(:,k)=2;
           T0_2(:,k)=-1;
           if mod(k,2)==0
               R1_2(:,k)=2;
               juzhen1=aa;
           else
               R1_1(:,k)=1;
               juzhen1=a;  
           end
       end
       
       if mod(k,2)==1
       barcer2=a*(X2(:,k)-X0(:,k));
       celiang2=a*(X2(:,k)-last_X2);
       else
       barcer2=aa*(X2(:,k)-X0(:,k));
       celiang2=aa*(X2(:,k)-last_X2);
       end
       
       s2bar(:,k)=g2*(barcer2);
       u2(:,k)=-K2*(barcer2)+u0(:,k)-miu*norm(barcer2)*s2bar(:,k)/(norm(s2bar(:,k))+0.0005);
       
       X2(:,k+1)=A*X2(:,k)+B*u2(:,k);       
       Eta2(:,1)=0.6;
       Eta2(:,k+1)=0.01*Eta2(:,k)+(Emax*barcer2'*phi*barcer2-(celiang2)'*phi0*(celiang2));

       if (celiang2)'*phi0*(celiang2)>=Emax*barcer2'*phi*barcer2+0*Eta2(:,k)
           Q2(:,k)=2;
           T2_1(:,k)=3;
           if (aa*(X2(:,k)-last_X2))'*(aa*(X2(:,k)-last_X2))>(a*(X2(:,k)-last_X2))'*(a*(X2(:,k)-last_X2))
               juzhen2=aa;
               R2_2(:,k)=2;
           else
               juzhen2=a;
               R2_1(:,k)=1;
           end
           last_X2=X2(:,k);
       elseif(celiang2)'*phi0*(celiang2)>=Emin*barcer2'*phi*barcer2
           last_X2=X2(:,k);
           Q2(:,k)=2;
           T2_2(:,k)=4;
           if mod(k,2)==1
               juzhen2=a;
               R2_1(:,k)=1;
           else
               juzhen2=aa;
               R2_2(:,k)=2;
           end
       end
       
       barcer3=juzhen2*(X3(:,k)-last_X2)+juzhen1*(X3(:,k)-last_X1);
       celiang3=juzhen2*(X3(:,k)-last_X3)+juzhen1*(X3(:,k)-last_X3);
       s3bar(:,k)=g3*(barcer3);
       u3(:,k)=-K3*(barcer3)+u0(:,k)-miu*norm(barcer3)*s3bar(:,k)/(norm(s3bar(:,k))+0.0005);

       X3(:,k+1)=A*X3(:,k)+B*u3(:,k);   
       Eta3(:,1)=1;
       Eta3(:,k+1)=0.01*(Eta3(:,k)+0.01*(Emax*barcer3'*phi*barcer3-(celiang3)'*phi0*(celiang3)));

       if (celiang3)'*phi0*(celiang3)>Emax*barcer3'*phi*barcer3+0*Eta3(:,k)
           Q3(:,k)=3;
           T3_1(:,k)=5;
           if (aa*(X3(:,k)-last_X3))'*(aa*(X3(:,k)-last_X3))>(a*(X3(:,k)-last_X3))'*(a*(X3(:,k)-last_X3))
               juzhen3=aa;
               R3_2(:,k)=2;
           else
               juzhen3=a;
               R3_1(:,k)=1;
           end
           last_X3=X3(:,k);
       elseif(celiang3)'*phi0*(celiang3)>Emin*barcer3'*phi*barcer3
           last_X3=X3(:,k);
           Q3(:,k)=3;
           T3_2(:,k)=6;
           if mod(k,2)==0
               juzhen3=aa;
               R3_2(:,k)=2;
           else
               juzhen3=a;
               R3_1(:,k)=1;
           end
       end

       barcer4=juzhen3*(X4(:,k)-last_X3);
       s4bar(:,k)=g4*(barcer4);
       u4(:,k)=-K4*(barcer4)+u0(:,k)-miu*norm(barcer4)*s4bar(:,k)/(norm(s4bar(:,k))+0.0005);

       X4(:,k+1)=A*X4(:,k)+B*u4(:,k);         
    end

    yizhixingwucha1=X0-X1;
    yizhixingwucha2=X0-X2;
    yizhixingwucha3=X1-X3+X2-X3;
    yizhixingwucha4=X3-X4;
    
    
figure(1)
plot3(2:N,X1(1,2:N),X1(2,2:N),'Color' ,' [0 0.4470 0.7410]','linewidth',4)
hold on
plot3(2:N,X2(1,2:N),X2(2,2:N),'Color' ,'[0.8500 0.3250 0.0980]','linewidth',4)
hold on
plot3(2:N,X3(1,2:N),X3(2,2:N),'Color' ,'[0.9290 0.6940 0.1250]','linewidth',4)
hold on
plot3(2:N,X4(1,2:N),X4(2,2:N),'Color' ,'[0.4660 0.6740 0.1880]','linewidth',4)
hold on
plot3(2:N,X0(1,2:N),X0(2,2:N),'Color' ,'[0.5 0.7 1]','linewidth',4)

xlabel('$Time(k)$','Interpreter','latex');
ylabel('$X_{i1}(k)$','Interpreter','latex')
zlabel('$X_{i2}(k)$','Interpreter','latex')
legend({'$agent 1$','$agent 2$','$agent 3$','$agent 4$','$leader$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',32);   


z0=zeros(N-1,1);
z1=ones(N-1,1);
z2=ones(N-1,1)*2;
z3=ones(N-1,1)*3;
z4=ones(N-1,1)*4;

figure(2)
plot3(2:N,z1,u1(2:N),'Color' ,' [0 0.4470 0.7410]','linewidth',1.2)
hold on
plot3(2:N,z2,u2(2:N),'Color' ,'[0.8500 0.3250 0.0980]','linewidth',1.2)
hold on
plot3(2:N,z3,u3(2:N),'Color' ,'[0.9290 0.6940 0.1250]','linewidth',1.2)
hold on
plot3(2:N,z4,u4(2:N),'Color' ,'[0.4660 0.6740 0.1880]','linewidth',1.2)
hold on
plot3(2:N,z0,u0(2:N),'Color' ,'[0.5 0.7 1]','linewidth',1.2)
xlabel('$Time(k)$','Interpreter','latex');
ylabel('Agent $i$','Interpreter','latex');
zlabel('$u_{i}(k)$','Interpreter','latex')
legend({'$agent 1$','$agent 2$','$agent 3$','$agent 4$','$leader$',},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',32);   

figure(3)
plot(2:N,s1bar(2:N),'Color' ,' [0 0.4470 0.7410]','linewidth',4)
hold on
plot(2:N,s2bar(2:N),'Color' ,'[0.8500 0.3250 0.0980]','linewidth',4)
hold on
plot(2:N,s3bar(2:N),'Color' ,'[0.9290 0.6940 0.1250]','linewidth',4)
hold on
plot(2:N,s4bar(2:N),'Color' ,'[0.4660 0.6740 0.1880]','linewidth',4)
xlabel('$Time(k)$','Interpreter','latex');
ylabel('$s_{i}(k)$','Interpreter','latex')
legend({'$\overline{s}_1$','$\overline{s}_2$','$\overline{s}_3$','$\overline{s}_4$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',32);   



figure(4)
plot(1:N,Q3(1:N),'*','Color' ,'[0.9290 0.6940 0.1250]','linewidth',1)
hold on
plot(1:N,Q2(1:N),'^','Color' ,'[0.8500 0.3250 0.0980]','linewidth',1)
hold on
plot(1:N,Q1(1:N),'p','Color' ,' [0 0.4470 0.7410]','linewidth',1)

xlabel('$Time(k)$','Interpreter','latex');
ylabel('$Tags$ $for$ $agents$','Interpreter','latex')
legend({'$agent 3$','$agent 2$','$agent 1$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',32);   
    

figure(5)
plot3( 2:N,z1,yizhixingwucha1(1,2:N),'Color' ,' [0 0.4470 0.7410]','linewidth',4)
hold on
plot3( 2:N,z2,yizhixingwucha2(1,2:N),'Color' ,'[0.8500 0.3250 0.0980]','linewidth',4)
hold on
plot3(2:N,z3, yizhixingwucha3(1,2:N),'Color' ,'[0.9290 0.6940 0.1250]','linewidth',4)
hold on
plot3(2:N,z4, yizhixingwucha4(1,2:N),'Color' ,'[0.4660 0.6740 0.1880]','linewidth',4)

xlabel('$Time(k)$','Interpreter','latex');
ylabel('Agent $i$','Interpreter','latex')
zlabel('$\bar{\xi}_{i}(k)$','Interpreter','latex')

legend({'$\bar{\xi}_1$','$\bar{\xi}_2$','$\bar{\xi}_3$','$\bar{\xi}_4$',},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',32); 


figure(6)
plot(1:N,T0_1(1:N),'*','Color' ,'[0 0 0]','linewidth',1)
hold on
plot(1:N,T0_2(1:N),'^','Color' ,'[0 0 0]','linewidth',1)
hold on
plot(1:N,T1_1(1:N),'*','Color' ,'r','linewidth',1)
hold on
plot(1:N,T1_2(1:N),'^','Color' ,'r','linewidth',1)
hold on
plot(1:N,T2_1(1:N),'*','Color' ,'g','linewidth',1)
hold on
plot(1:N,T2_2(1:N),'^','Color' ,'g','linewidth',1)
hold on
plot(1:N,T3_1(1:N),'*','Color' ,'b','linewidth',1)
hold on
plot(1:N,T3_2(1:N),'^','Color' ,'b','linewidth',1)
xlabel('$Time(k)$','Interpreter','latex');
legend({'$WTOD$','$RR$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',38);

figure(7)
plot(1:N,T2_1(1:N),'*','Color' ,'[0.48 0.67 0.87]','linewidth',1)
hold on
plot(1:N,T2_2(1:N),'^','Color' ,'[0.58 0.77 0.63]','linewidth',1)

xlabel('$Time(k)$','Interpreter','latex');
ylabel('$Agent2$ $scheduling$ $protocols$','Interpreter','latex')
legend({'$WTOD$','$RR$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',38);

figure(8)
plot(1:N,T3_1(1:N),'*','Color' ,'[0.48 0.67 0.87]','linewidth',1)
hold on
plot(1:N,T3_2(1:N),'^','Color' ,'[0.58 0.77 0.63]','linewidth',1)

xlabel('$Time(k)$','Interpreter','latex');
ylabel('$Agent3$ $scheduling$ $protocols$','Interpreter','latex')
legend({'$WTOD$','$RR$'},'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',38);







