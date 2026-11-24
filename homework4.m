syms L theta1 theta2 theta3
assume(theta1>=-pi & theta1<=pi & theta2>=-pi & theta2<=pi & L>0) 
screw1 = [0 0 1 0 -L 0];
screw2 = [0 1 0 L 0 0];
screw3 = [0 0 0 1 0 0];
s1 = bracket(screw1);
s2 = bracket(screw2);
s3 = bracket(screw3);
M = [1 0 0 2*L; 0 1 0 L; 0 0 1 -L; 0 0 0 1];
Tee = expm(s1*theta1)* expm(s2*theta2)* expm(s3*theta3)*M;
Teesimple=simplify(Tee)
pee = Tee(1:3, 4);
peeSimple=simplify(pee)
eqn1  = Tee == [0 1 0 2*L; 0 0 -1 0;-1 0 0 -3*L; 0 0 0 1];
vars = [theta1, theta2, theta3];
sol1 = solve(eqn1, vars, "Real",true, "ReturnConditions",true, PrincipalValue=true)

%homogenous Transform 
c1=cos(theta1); s1=sin(theta1);
c2=cos(theta2); s2=sin(theta2);

ts1  = [c1 -s1 0 L;
        s1 c1 0 0;
        0 0 1 0;
        0 0 0 1];   % Trans(L,0,0)*Rz(θ1)
t12  = [c2 0 s2 0;
        0 1 0 L;
        -s2 0 c2 -L;
        0 0 0 1]; % Trans(0,L,-L)*Ry(θ2)
t23  = [eye(3) [L;0;0]; 0 0 0 1];                      % Trans(L,0,0) to joint 3
t3ee = [eye(3) [theta3;0;0]; 0 0 0 1];                 % prismatic slide along x

Tsee_check = simplify(ts1*t12*t23*t3ee)

eqn  = Tsee_check == [0 1 0 2*L; 0 0 -1 0;-1 0 0 -3*L; 0 0 0 1];
% Tee = 
vars = [theta1, theta2, theta3];

sol = solve(eqn, vars, "Real",true, "ReturnConditions",true, PrincipalValue=true)


function S_matrix = bracket(S)
    S_matrix = [ 0 -S(3) S(2) S(4);
                S(3) 0 -S(1) S(5);
                -S(2) S(1) 0 S(6);
                0 0 0 0];
end

I = eye(3);
T01 = expm(s1*theta1)* [I, [L;0;0] ; 0 0 0 1];
T02 = expm(s1*theta1)* expm(s2*theta2) * [I, [L;0;-L] ; 0 0 0 1];
T03 = expm(s1*theta1)* expm(s2*theta2) * expm(s3*theta3)* [I, [0;L;-L] ; 0 0 0 1];


% te =[theta1;theta2;theta3];
% Jv = jacobian(pee,te);
% size(Jv);
% Jvsimple=simplify(Jv)
% Jw = [T01(1:3,3), T02(1:3, 2), zeros(3,1)];
% Jwsimple=simplify(Jw)
% latJv = latex(Jvsimple);
% latJw = latex(Jwsimple);
% latPee = latex(peeSimple);

A = [c1 s1 0; -s1 c1 0; 0 0 1];
B = [c2 s2 L; -s2 c2 0; 0 0 1];
C = [1 0 L; 0 1 0; 0 0 1];
Trr = A*B*C;
Prr = Trr(1:2,3);
Jv = jacobian(simplify(Prr),[theta1,theta2]);
detJv= det(Jv)

eqn  = detJv == 0;
vars = [theta1, theta2];
solSing = solve(eqn, vars, "Real",true, "ReturnConditions",true)
sqrt(det(Jv*Jv'))
