clear; clc; close all;
syms m1 m2 r1 r2 l1 l2 I1 I2 q1 q2 dq1 dq2 ddq1 ddq2 u1 u2 g 'real'
syms t 
m1=0.75; m2=0.75; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.063; I1= 0.063;

% Generation of a cubic trajectory for the 1st joint
tspan = [0, 10];
x0 = [pi; 0];
xf = [0; 0];
a1 = CubicTraj(tspan, x0, xf);

% Generation of a cubic trajectory for the 2nd joint
tspan = [0, 10];
x0 = [pi/2; 0];
xf = [0; 0];
a2 = CubicTraj(tspan, x0, xf);

% Eigenvalue placement method
p = [-3, -3, -4, -4];
A = [0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0];
B = [0, 0; 0, 0; 1, 0; 0, 1];
K = place(A, B, p);

% Finding Lyuapunov P matrix
Acl = [0, 0, 1, 0; 0, 0, 0, 1; -12, 0, -7, 0; 0, -12, 0, -7];
Q = eye(4);
P = lyap(Acl', Q);

% Boundary layer
phi = 0.006;

% Simulation of the system 
x0 = [deg2rad(200),deg2rad(125),0, 0];
[T,X] = ode45(@(t,x) ode2linkTracking(t,x, K, P, phi), tspan, x0);
% [T,X] = ode45(@ode2linkTracking, tspan, x0);

q1_d = pi - (3*pi.*T.^2)/100 + (pi.*T.^3)/500;
q1dot_d = - (6*pi.*T)/100 + (3.*T.^2*pi)/500;
q1ddot_d = - (6*pi)/100 + (6.*T*pi)/500;

q2_d =  pi/2 - (3*pi.*T.^2)/200 + (pi.*T.^3)/1000;
q2dot_d = - (6*pi.*T)/200 + (3.*T.^2*pi)/1000;
q2ddot_d =- (6*pi)/200 + (6.*T*pi)/1000;

% Store the torque values
Tau = [];
for index = 1:length(T)
   time = T(index);
   x = X(index, :).';
   [~,tau] = ode2linkTracking(time, x, K, P, phi);
   Tau = [Tau, tau];
end

figure(1)
subplot(3,2,1);
plot(T, q1_d);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(T,X(:,1),'r');
xlabel('t', 'FontSize',14)
ylabel('theta1','FontSize',14);
legend('desired','actual')
hold off

subplot(3,2,2);
plot(T, q1dot_d);
xlabel('t', 'FontSize',14)
ylabel('q1_d','FontSize',14);
hold on
plot(T,X(:,3),'r');
xlabel('t', 'FontSize',14)
ylabel('theta1 dot','FontSize',14)
legend('desired','actual')
hold off

subplot(3,2,3);
plot(T, q2_d);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(T,X(:,2),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2','FontSize',14)
legend('desired','actual')
hold off

subplot(3,2,4);
plot(T, q2dot_d);
xlabel('t', 'FontSize',14)
ylabel('q2_d','FontSize',14);
hold on
plot(T,X(:,4),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2 dot','FontSize',14)
legend('desired','actual')
hold off

subplot(3, 2, 5);
plot(T, Tau(1, :),'b');
xlabel('t', 'FontSize', 14)
ylabel('u1', 'FontSize', 14);

subplot(3, 2, 6);
plot(T, Tau(2, :),'b');
xlabel('t', 'FontSize', 14)
ylabel('u2', 'FontSize', 14);
