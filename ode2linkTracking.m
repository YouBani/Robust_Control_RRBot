function [dX, tau] = ode2linkTracking(t,x)
m1_nom=0.75; m2_nom=0.75; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I1_nom= 0.063; I2_nom= 0.063; m1=1; m2=1;I2= 0.084; I1= 0.084;

dX= zeros(4,1);
x=num2cell(x);
[theta1, theta2, theta1_dot, theta2_dot] = deal(x{:});

a = I1_nom + I2_nom + m1_nom*r1^2 + m2_nom*(l1^ 2 + r2^2);
b = m2_nom*l1*r2;
d = I2_nom + m2_nom*r2^2;

Mmat= [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
Cmat= [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
Gmat= [-m1_nom*g*r1*sin(theta1)-m2_nom*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2_nom*g*r2*sin(theta1+theta2)];
% Mmat= [a+2*b*cos(x(2)), d+b*cos(x(2)); d+b*cos(x(2)), d];
% Cmat= [-b*sin(x(2))*x(4), -b*sin(x(2))*(x(3)+x(4)); b*sin(x(2))*x(3),0];
% Gmat= [-m1*g*r1*sin(x(1))-m2*g*(l1*sin(x(1))+r2*sin(x(1)+x(2))); -m2*g*r2*sin(x(1)+x(2))];
% Desired trajectories
q1_d = pi - (3*pi*t^2)/100 + (pi*t^3)/500;
q1dot_d = - (6*pi*t)/100 + (3*t^2*pi)/500;
q1ddot_d = - (6*pi)/100 + (6*t*pi)/500;

q2_d =  pi/2 - (3*pi*t^2)/200 + (pi*t^3)/1000;
q2dot_d = - (6*pi*t)/200 + (3*t^2*pi)/1000;
q2ddot_d = - (6*pi)/200 + (6*t*pi)/1000;

% Feedback linearization control 
K = [12, 0, 7, 0; 0, 12, 0, 7];
B = [0, 0; 0, 0; 1, 0; 0, 1];
P = [1.2202, 0, 0.0417, 0; 0, 1.2202, 0, 0.0417; 0.0417, 0, 0.0774, 0; 0, 0.0417, 0, 0.0774];
% P = [3.631, 3.631, 0.4167, 0.4167; 0, 0, 0, 0; 0.4167, 0.4167, 0.0595, 0.0595; 0, 0, 0, 0];

x = [theta1 - q1_d; theta2 - q2_d; theta1_dot - q1dot_d; theta2_dot - q2dot_d];
rau = 15;
phi = 0.006;

if norm(B'*P*x) > phi
    vr = -rau * (B'*P*x) / norm(B'*P*x);
else
    vr = -rau * (B'*P*x) / phi;
end

% if norm(B'*P*x) > 0
%    vr = -rau * (B'*P*x) / norm(B'*P*x);
% else
%     vr = 0;
% end

% v = - K*x + [q1ddot_d; q2ddot_d] + vr;
v = - K*x + [q1ddot_d; q2ddot_d];
tau = Mmat*v + Cmat*[theta1_dot; theta2_dot] + Gmat;

u1 = [tau(1,1)];
u2 = [tau(2,1)];

dX(1) = theta1_dot;

dX(2) = theta2_dot;

% ddq = Mmat\(tau - Cmat*[theta1_dot; theta2_dot] - Gmat);
% dX(3) = ddq(1, 1);
% dX(4) = ddq(2, 1);

dX(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dX(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

end