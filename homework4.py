import sympy as sp

# Define symbolic variables
L, theta1, theta2, theta3 = sp.symbols('L theta1 theta2 theta3', real=True)

# 4x4 homogeneous transformation matrices
Ts1 = sp.Matrix([
    [1, 0, 0, L],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

T12 = sp.Matrix([
    [sp.cos(theta1), -sp.sin(theta1), 0, 0],
    [sp.sin(theta1),  sp.cos(theta1), 0, L],
    [0, 0, 1, -L],
    [0, 0, 0, 1]
])

T23 = sp.Matrix([
    [ sp.cos(theta2), 0, sp.sin(theta2), L],
    [0, 1, 0, 0],
    [-sp.sin(theta2), 0, sp.cos(theta2), 0],
    [0, 0, 0, 1]
])

T3ee = sp.Matrix([
    [1, 0, 0, theta3],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
])

# Homogeneous transformation multiplication
Tsee = Ts1 * T12 * T23 * T3ee

print(sp.pretty(Tsee))
print("Simplified:")
print(sp.pretty(Tsee.applyfunc(sp.trigsimp)))

Tcurrent = sp.Matrix([
    [0, 1, 0, 2*L],
    [0, 0, -1, 0],
    [-1, 0, 0, -3*L],
    [0, 0, 0, 1]
])

# Solve Tsee = Tcurrent for theta1, theta2, and theta3
solution = sp.solve(
    [Tsee[i, j] - Tcurrent[i, j]
     for i in range(4)
     for j in range(4)],
    (theta1, theta2, theta3),
    dict=True
)

print("Tcurrent:")
print(sp.pretty(Tcurrent))
print("Solutions:")
for joint_values in solution:
    print(joint_values)
    print(len(joint_values))
    
