syms L theta1 theta2 theta3
assume(theta1>=-pi & theta1<=pi & theta2>=pi & theta2<=pi & L>0) 
screw1 = [0 0 1 0 -L 0];
screw2 = [0 1 0 L 0 0];
screw3 = [0 0 0 1 0 0];
s1 = bracket(screw1);
s2 = bracket(screw2);
s3 = bracket(screw3);

M = [1 0 0 2*L; 0 1 0 L; 0 0 1 -L; 0 0 0 1];
Tee = expm(s1*theta1)* expm(s2*theta2)* expm(s3*theta3)*M;
Teesimple=simplify(Tee)
I = eye(3);
T01 = expm(s1*theta1)* [I, [L;0;0] ; 0 0 0 1];
T02 = expm(s1*theta1)* expm(s2*theta2) * [I, [L;0;-L] ; 0 0 0 1];
T03 = expm(s1*theta1)* expm(s2*theta2) * expm(s3*theta3)* [I, [0;L;-L] ; 0 0 0 1];

pee = Tee(1:3, 4);
peeSimple=simplify(pee)

eqn  = pee == [2*L; 0;-3*L];
% Tee = 
vars = [theta1, theta2, theta3];

sol = solve(eqn, vars, "Real",true, "ReturnConditions",true, PrincipalValue=true)

% te =[theta1;theta2;theta3];
% Jv = jacobian(pee,te);
% size(Jv);
% Jvsimple=simplify(Jv)
% Jw = [T01(1:3,3), T02(1:3, 2), zeros(3,1)];
% Jwsimple=simplify(Jw)
% latJv = latex(Jvsimple);
% latJw = latex(Jwsimple);
% latPee = latex(peeSimple);

ts1= [eye(4)+[[0,0,0,L];zeros(3,4)]];
t12= rot

function S_matrix = bracket(S)
    S_matrix = [ 0 -S(3) S(2) S(4);
                S(3) 0 -S(1) S(5);
                -S(2) S(1) 0 S(6);
                0 0 0 0];
    end