%general params
J=1;
K=1;
n=1;
tau=1;
%open ("motor_model.slx");
v=1;
g=9.81;
kT = 0.17;
feff=2.4*10^-5;
ni = 53;
Jm = 1320*10^-7;

%general robot dims
L1 = 0.67;
r1 = 0.04;
L2 = 1.7;
b2 = 0.22;
L3 = 1.65;
b3 = 0.18;
a4 = 0.98;
m1 = 4.9;
m2 = 8.1;
m3 = 4.9;
m4 = 2.2;
dd2 = 0.34;

%joint params
q1 = 0;
q2 = (pi/6+3*pi/4)/2;
q3 = 2.175;
q4 = 7/8*pi;
q = [q1,q2,q3,q4];

%Interatia Matrices
D1 = [1/12*m1*(3*r1^2+L1^2) 0 0;
    0 0.5*m1*r1^2 0;
    0 0 1/12*m1*(3*r1^2+L1^2)];
D2 = [1/12*m2*(b2^2+L2^2) 0 0;
    0 1/12*m2*(b2^2+L2^2) 0;
    0 0 1/12*m2*(b2^2+b2^2)];
D3 = [(1/12)*m3*(b3^2+L3^2),0,0;
    0,(1/12)*m3*(b3^2+b3^2),0;
    0,0,(1/12)*m3*(b3^2+L3^2)];
D4 = [0,0,0;
    0,(1/12)*m4*a4^2,0;
    0,0,(1/12)*m4*a4^2];

%Intertia Matrix indices
I4 = D4(2,2);
I2 = D2(1,1);
I2zz = D2(3,3);
I3 = D3(1,1);
I3yy = D3(2,2);
I1yy = D1(2,2);

%Dynamics calculations
K0 = m2*(0.5*L2-dd2)^2+m3*(0.5*L3)^2;
K1 = 0.5*(I4+0.25*m4*a4^2);
K2 = 0.5*(I2-I2zz+I3-I3yy+K0);
K3 = I2+I3+K0+2*K1;
K4 = m3+m4;
K5 = 0.5*(2*I1yy+I2zz+I3yy+K3);
    
%lower and upper bounds for fmin
lb = [-pi,pi/6,1.35,-pi/2];
ub = [pi,3*pi/4,3,5*pi/4];

%functions to calculate J_eff
f_11 = @(q)-(K5+(0.5*K4*q(3)^2-0.5*m3*L3*q(3))-(K2+(0.5*K4*q(3)^2-0.5*m3*L3*q(3)))*cos(2*q(2))+K1*cos(2*(q(2)+q(4)))-(m4*a4*q(3))*cos(q(2)+q(4))*sin(q(2)));
f_22 = @(q)-(K3 + 2*(0.5*K4*q(3)^2-0.5*m3*L3*q(3))+(m4*a4*q(3))*sin(q(4)));
%%
%problem 12

x_11 = fmincon(f_11,q,[],[],[],[],lb,ub);
x_22 = fmincon(f_22,q,[],[],[],[],lb,ub);

w_n = 15;
Zeta = 1;
q_r = 0.35;
f = @(i)(1/ni^2)*i+Jm;

J_eff = [f(-f_11(x_11)), f(-f_22(x_22)),f(K4),f(2*K1)];

kd_f = @(i)2*Zeta/w_n*K_p(i)-feff/kT;
K_p = [w_n^2*J_eff(1)/kT,w_n^2*J_eff(2)/kT,w_n^2*J_eff(3)/kT,w_n^2*J_eff(4)/kT];
K_d = [kd_f(1),kd_f(2),kd_f(3),kd_f(4)];
%% Problem 13
h_q = [0 0.5*m4*a4*cos(q2+q4)+(m2*(dd2-0.5*L2)+m3*(0.5*L3-q3)-m4*q3)*sin(q2) (m3+m4)*cos(q2) 0.5*a4*m4*cos(q2+q4)];

fhq_2 = @(q)-(0.5*m4*a4*cos(q(2)+q(4))+(m2*(dd2-0.5*L2)+m3*(0.5*L3-q(3))-m4*q(3))*sin(q(2)));
fhq_3 = @(q)-((m3+m4)*cos(q(2)));
fhq_4 = @(q)-(0.5*a4*m4*cos(q(2)+q(4)));

h_2 = fmincon(fhq_2,q,[],[],[],[],lb,ub);
h_3 = fmincon(fhq_3,q,[],[],[],[],lb,ub);
h_4 = fmincon(fhq_4,q,[],[],[],[],lb,ub);

q_max = [max([h_2;h_3;h_4],[],1)];
Tl = [0 fhq_2(q_max) fhq_3(q_max) fhq_4(q_max)];
%%


K_I = [w_n^2/kT,w_n^2/kT,w_n^2/kT,w_n^2/kT];
K_P = [(2*Zeta*w_n+w_n^2*abs(J_eff(1)))/kT,(2*Zeta*w_n+w_n^2*abs(J_eff(2)))/kT,(2*Zeta*w_n+w_n^2*abs(J_eff(3)))/kT,(2*Zeta*w_n+w_n^2*abs(J_eff(4)))/kT];
K_D = [(2*Zeta*w_n*abs(J_eff(1))+1-feff)/kT,(2*Zeta*w_n*abs(J_eff(2))+1-feff)/kT,(2*Zeta*w_n*abs(J_eff(3))+1-feff)/kT,(2*Zeta*w_n*abs(J_eff(4))+1-feff)/kT];
%%
