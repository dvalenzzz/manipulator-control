syms L theta1 theta2 theta3 theta4
assume(theta1>=-2*pi & theta1<=2*pi & theta2>=-2*pi & theta2<=2*pi & theta3>=-2*pi & theta3<=2*pi & theta4>=-2*pi & theta4<=2*pi)
screw1 = [0 0 -1 0 0 0];
screw2 = [0 -1 0 2*L 0 -L];
screw3 = [-1 0 0 0 -L 0];
screw4 = [0 0 0 1 0 0];
s1 = bracket(screw1);
s2 = bracket(screw2);
s3 = bracket(screw3);
s4 = bracket(screw4);

M = [0 0 1 2*L; 0 1 0 0; -1 0 0 L; 0 0 0 1];
Tee = expm(s1*theta1)* expm(s2*theta2)* expm(s3*theta3)* expm(s4*theta4)*M;
Teesimple=simplify(Tee)
I = eye(3);
T01 = expm(s1*theta1)* [I, [0;0;2*L] ; 0 0 0 1];
T02 = expm(s1*theta1)* expm(s2*theta2) * [I, [L;0;2*L] ; 0 0 0 1];
T03 = expm(s1*theta1)* expm(s2*theta2) * expm(s3*theta3)* [I, [L;0;L] ; 0 0 0 1];

pee = Tee(1:3, 4);
peeSimple=simplify(pee)
te =[theta1;theta2;theta3;theta4];
Jv = jacobian(pee,te);
size(Jv);
Jvsimple=simplify(Jv)
Jw = [T01(1:3,3), T02(1:3, 2), T03(1:3,1), [0;0;0]];
Jwsimple=simplify(Jw)
latJv = latex(Jvsimple);
latJw = latex(Jwsimple);
latPee = latex(peeSimple);

function S_matrix = bracket(S)
    S_matrix = [ 0 -S(3) S(2) S(4);
                S(3) 0 -S(1) S(5);
                -S(2) S(1) 0 S(6);
                0 0 0 0];
    end