clc;

clear all;

close all;

T=45;% Simulation Time

K=5;

P=10;

Ki=0.005;

Sum=0;

h=0.01;

f=0.05;

I=diag([100;75;80]);

sigma_BN=[0.1;0.2;-0.1];

omega_BN=(pi/180)*[3;1;-2];

x=[sigma_BN;omega_BN];

SigArr=[];

for t=0:h:T

      sigma_BN=x(1:3);

      omega_BN=x(4:6);

      sigma_RN=[0.2*sin(f*t);0.3*cos(f*t);-0.3*sin(f*t);];

      sigma_BR=addMRP(-sigma_RN,sigma_BN);

      if norm(sigma_BR)>1

          sigma_BR=-sigma_BR/norm(sigma_BR)^2;

      end

      SigArr=[SigArr sigma_BR];

      C=MRP2C(sigma_BR);

      sigma_RN_dot=[0.2*f*cos(f*t);-0.3*f*sin(f*t);-0.3*f*cos(f*t);];

      omega_RN=4/(1+norm(sigma_RN)^2)^2*BmatMRP(sigma_RN)'*sigma_RN_dot;

      if t==0

         delta_omega_0=omega_BN-C*omega_RN;

      end

      %The following terms is the analytical solution of omega_RN_dot

      omega_RN_dot=[(6*f*sin(f*t)*((3*f*cos(f*t))/5 - (3*f*cos(f*t)^2)/25 + (3*f*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f^2*cos(f*t)*((3*sin(f*t))/5 - (3*cos(f*t)*sin(f*t))/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (4*f*cos(f*t)*((4*f*cos(f*t)*sin(f*t))/25 + (9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (4*f^2*sin(f*t)*((2*sin(f*t)^2)/25 - (9*abs(cos(f*t))^2)/100 - (13*abs(sin(f*t))^2)/100 + 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (6*f*cos(f*t)*((3*f*sin(f*t))/5 - (6*f*cos(f*t)*sin(f*t))/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (6*f^2*sin(f*t)*((3*cos(f*t))/5 + (3*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (12*f*sin(f*t)*((3*sin(f*t))/5 - (3*cos(f*t)*sin(f*t))/25)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) + (8*f*cos(f*t)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50)*((2*sin(f*t)^2)/25 - (9*abs(cos(f*t))^2)/100 - (13*abs(sin(f*t))^2)/100 + 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) + (12*f*cos(f*t)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50)*((3*cos(f*t))/5 + (3*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3);

                    (4*f*cos(f*t)*((3*f*cos(f*t))/5 + (3*f*cos(f*t)^2)/25 - (3*f*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (6*f*cos(f*t)*((2*f*cos(f*t))/5 - (9*f*cos(f*t)^2)/50 + (9*f*sin(f*t)^2)/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (4*f^2*sin(f*t)*((3*sin(f*t))/5 + (3*cos(f*t)*sin(f*t))/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f^2*sin(f*t)*((2*sin(f*t))/5 - (9*cos(f*t)*sin(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f*sin(f*t)*((9*f*cos(f*t)*sin(f*t))/25 - (9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 + (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f^2*cos(f*t)*((9*abs(cos(f*t))^2)/100 + (13*abs(sin(f*t))^2)/100 - (9*cos(f*t)^2)/50 - 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (8*f*cos(f*t)*((3*sin(f*t))/5 + (3*cos(f*t)*sin(f*t))/25)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) - (12*f*cos(f*t)*((2*sin(f*t))/5 - (9*cos(f*t)*sin(f*t))/50)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) + (12*f*sin(f*t)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50)*((9*abs(cos(f*t))^2)/100 + (13*abs(sin(f*t))^2)/100 - (9*cos(f*t)^2)/50 - 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3);

                    (6*f*sin(f*t)*((2*f*cos(f*t))/5 + (9*f*cos(f*t)^2)/50 - (9*f*sin(f*t)^2)/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f^2*cos(f*t)*((2*sin(f*t))/5 + (9*cos(f*t)*sin(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (6*f*cos(f*t)*((9*f*cos(f*t)*sin(f*t))/25 + (9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (6*f^2*sin(f*t)*((9*sin(f*t)^2)/50 - (9*abs(cos(f*t))^2)/100 - (13*abs(sin(f*t))^2)/100 + 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (4*f*cos(f*t)*((3*f*sin(f*t))/5 + (6*f*cos(f*t)*sin(f*t))/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) - (4*f^2*sin(f*t)*((3*cos(f*t))/5 - (3*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^2) + (12*f*sin(f*t)*((2*sin(f*t))/5 + (9*cos(f*t)*sin(f*t))/50)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) - (12*f*cos(f*t)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50)*((9*sin(f*t)^2)/50 - (9*abs(cos(f*t))^2)/100 - (13*abs(sin(f*t))^2)/100 + 1))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3) + (8*f*cos(f*t)*((9*f*sign(cos(f*t))*sin(f*t)*abs(cos(f*t)))/50 - (13*f*abs(sin(f*t))*sign(sin(f*t))*cos(f*t))/50)*((3*cos(f*t))/5 - (3*sin(f*t)^2)/25))/(5*((13*abs(sin(f*t))^2)/100 + (9*abs(cos(f*t))^2)/100 + 1)^3);];

      u=-K*sigma_BR-(P*eye(3)+P*Ki*I)*(omega_BN-C*omega_RN)-K*P*Ki*Sum+P*Ki*I*delta_omega_0+I*(C*omega_RN_dot-skew(omega_BN)*C*omega_RN)+skew(omega_BN)*I*omega_BN;

      K1=dotx(x,u);

      K2=dotx(x+0.5*h*K1,u);

      K3=dotx(x+0.5*h*K2,u);

      K4=dotx(x+h*K3,u);

      x=x+h/6*(K1+2*K2+2*K3+K4);

      %x=x+dotx(x,u)*h;

      if norm(x(1:3))>1

          x(1:3)=-x(1:3)/norm(x(1:3))^2;

      end

      Sum=Sum+sigma_BR*h;% integral temrs computation

end

delta=norm(sigma_BR)


%plot figure

i=0:h:T;

plot(i,SigArr(1,:),'r',i,SigArr(2,:),'g',i,SigArr(3,:),'b');

title('\sigmabr')

grid


function f=skew(x)


f=[0,-x(3),x(2);

     x(3),0,-x(1);

     -x(2),x(1),0;];


end


function y=dotx(x,u)

sigma=x(1:3);

omega=x(4:6);

L=[0.5;-0.3;0.2];

B=(1-norm(sigma)^2)*eye(3)+2*skew(sigma)+2*sigma*sigma';

sigmadot=0.25*B*omega;

I=diag([100;75;80]);

omegadot=inv(I)*(-skew(omega)*I*omega+u+L);


y=[sigmadot;omegadot];

end