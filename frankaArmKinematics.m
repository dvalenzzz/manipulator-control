% syms theta1 theta2 theta3 theta4 theta5 theta6 theta7
% assume(theta1>=-2*pi & theta1<=2*pi & theta2>=-2*pi & theta2<=2*pi & theta3>=-2*pi & theta3<=2*pi & theta4>=-2*pi & theta4<=2*pi & theta5>=-2*pi & theta5<=2*pi & theta6>=-2*pi & theta6<=2*pi & theta7>=-2*pi & theta7<=2*pi)
if ispc
    % Windows path (Update this folder name to your Windows CasADi download)
    addpath(genpath('D:\casadi-3.7.2-windows64-matlab2018b'));
elseif isunix
    addpath(genpath('/home/dario/dualArm/casadi-3.7.2-linux64-matlab2018b'));  % or wherever you installed it
end
import casadi.*
theta = casadi.SX.sym('theta', 7);
theta1 = theta(1);
theta2 = theta(2);
theta3 = theta(3);
theta4 = theta(4);
theta5 = theta(5);
theta6 = theta(6);
theta7 = theta(7);
T01= homoTransform(0,0,0.333,0);%assumes base flange is the {s} frame origin
T12= homoTransform(0,0,0,-pi/2); %information gathered from https://github.com/frankarobotics/franka_description/blob/main/robots/fr3v2/kinematics.yaml
T23= homoTransform(0,-0.316,0,pi/2);
T34= homoTransform(0.0825,0,0,pi/2);
T45= homoTransform(-0.0825,0.384,0,-pi/2);
T56= homoTransform(0,0,0,pi/2);
T67= homoTransform(0.088,0,0,pi/2);
T7ee= homoTransform(0,0,0.107,0);
incrementalT = {T01;T12;T23;T34;T45;T56;T67;T7ee};
n = numel(incrementalT);
M = eye(4);
T0i= cell(1,n);
Mstore = cell(1,n);
for k = 1:n
M = M * incrementalT{k};
T0i{k} = M;
Mstore{k}=M;
end
% Mlocal=eye(4);
% Mstore = cell(1,n);
% for k = 1:n
%     Mstore{k} = MlocalincrementalT{k};
% end
Melbow = Mstore{4};
% robot = importrobot("fr3v2.urdf");
%axis of rotation must be in stationary frame which is found by
%RsiLocalAxis@jointi
w = cell(1,7);
q = cell(1,7);
s = cell(1,7);
S_matrix = cell(1,7);
for k = 1:7
w{k} = T0i{k}(1:3,1:3)*[0;0;1]; %local axis of rotation.
q{k} = (eye(3)-w{k}*w{k}')*T0i{k}(1:3,4); %minimum distance from the {s} origin to the k-joint
s{k} = [w{k};-cross(w{k},q{k})];
S_matrix{k} = bracket(s{k});
end
%now that we have the screw matrices and zero configuration matrix M we can do product of exponentials
%making sure to include the correct sign for axis of rotation.
s1 = S_matrix{1};s2 = S_matrix{2};s3 = S_matrix{3};s4 = S_matrix{4};s5 = S_matrix{5};s6 = S_matrix{6};s7 = S_matrix{7};
Telbow = expmCAS(s1,theta1) *expmCAS(s2,theta2)* expmCAS(s3,theta3)*expmCAS(s4,theta4)*Melbow; %TODO: make sure doing this at the correct joint. Prob. just assume a cylinder centered at that point the size of the elbow.

T01sym = Mstore{1}; % no joints before frame 1, so it's constant
T02sym = expmCAS(s1,theta1)*Mstore{2};
T03sym = expmCAS(s1,theta1)*expmCAS(s2,theta2)*Mstore{3};
T04sym = expmCAS(s1,theta1)*expmCAS(s2,theta2)*expmCAS(s3,theta3)*Mstore{4};
T05sym = expmCAS(s1,theta1)*expmCAS(s2,theta2)*expmCAS(s3,theta3)*expmCAS(s4,theta4)*Mstore{5};
T06sym = expmCAS(s1,theta1)*expmCAS(s2,theta2)*expmCAS(s3,theta3)*expmCAS(s4,theta4)*expmCAS(s5,theta5)*Mstore{6}; %you can also replace *theta1-theta7 with casadi symbols.
T07sym = expmCAS(s1,theta1) *expmCAS(s2,theta2) *expmCAS(s3,theta3)*expmCAS(s4,theta4)*expmCAS(s5,theta5)*expmCAS(s6,theta6)*Mstore{7}; %you can also replace theta1-theta7 with casadi symbols.
Tee = expmCAS(s1,theta1) *expmCAS(s2,theta2) *expmCAS(s3,theta3)*expmCAS(s4,theta4)*expmCAS(s5,theta5)*expmCAS(s6,theta6)*expmCAS(s7,theta7)*Mstore{8}; %you can also replace theta1-theta7 with casadi symbols.
T0isym = {T01sym, T02sym, T03sym, T04sym, T05sym, T06sym, T07sym, Tee};

pee = Tee(1:3, 4);
Jv = jacobian(pee,theta);
size(Jv);
Jw = [T0isym{1}(1:3,3), T0isym{2}(1:3, 3), T0isym{3}(1:3,3), T0isym{4}(1:3,3), T0isym{5}(1:3,3),T0isym{6}(1:3,3),T0isym{7}(1:3,3)];
size(Jw);
J = [Jw;Jv];
fkFun = casadi.Function('fkFun', {theta}, {Tee});
JFun  = casadi.Function('JFun',  {theta}, {J});

%AI smoke test
%% Smoke test + comparison (double-valued wrappers around the CasADi Functions)
fk = @(x) full(fkFun(x));   % DM -> double
Jf = @(x) full(JFun(x));

% finite-difference check of the Jacobian
h = 1e-6; q = [0 -pi/4 0 -3*pi/4 0 pi/2 pi/4]' + 0.3*randn(7,1);
T = fk(q); Jfd = zeros(6,7);
for i = 1:7
    dq = zeros(7,1); dq(i) = h;
    Tp = fk(q+dq); Tm = fk(q-dq);
    Vs = (Tp - Tm)/(2*h) / T;
    Jfd(:,i) = [Vs(3,2); Vs(1,3); Vs(2,1); (Tp(1:3,4) - Tm(1:3,4))/(2*h)];
end
fprintf("J     vs FD   : %.2e\n", norm(Jf(q) - Jfd));

%% Compare against the Robotics System Toolbox Panda model
panda  = loadrobot("frankaEmikaPanda", DataFormat="column");
eeName = "panda_link8";   % flange, same as Mstore{8}
qReady = [0 -pi/4 0 -3*pi/4 0 pi/2 pi/4]';

% 1) zero configuration: M
Tp0 = getTransform(panda, zeros(9,1), eeName);
fprintf("M     vs Panda: %.2e\n", norm(Mstore{8} - Tp0));

% 2) FK and Jacobian at several random poses
nTests = 20; errFK = zeros(1,nTests); errJ = zeros(1,nTests);
for t = 1:nTests
    qt = qReady + 0.5*randn(7,1);
    qPanda = [qt; 0; 0];   % Panda model has 2 finger joints
    errFK(t) = norm(fk(qt) - getTransform(panda, qPanda, eeName));
    Jp = geometricJacobian(panda, qPanda, eeName);   % [w; pdot], base frame
    errJ(t) = norm(Jf(qt) - Jp(:,1:7));
end
fprintf("FK    vs Panda: max %.2e over %d poses\n", max(errFK), nTests);
fprintf("J     vs Panda: max %.2e over %d poses\n", max(errJ), nTests);

function sequentialTransform = homoTransform(x,y,z,r) %currently only valid for rotations about the local x-axis. Full support should include product of RxRyRz

sequentialTransform = [1 0 0 x;
0 cos(r) -sin(r) y;
0 sin(r) cos(r) z;
0 0 0 1];

end
function S_matrix = bracket(S)
S_matrix = [ 0 -S(3) S(2) S(4);
S(3) 0 -S(1) S(5);
-S(2) S(1) 0 S(6);
0 0 0 0];
end
function T = expmCAS(Smat, th) 
% exp(Smat*th) for a revolute screw (unit w) - closed form, works with CasADi SX
W = Smat(1:3,1:3);   % [w]
v = Smat(1:3,4);
R = eye(3) + sin(th)*W + (1-cos(th))*(W*W);
G = eye(3)*th + (1-cos(th))*W + (th-sin(th))*(W*W);
T = [R, G*v; 0 0 0 1];
end