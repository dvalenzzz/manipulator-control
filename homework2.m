Wx = [0 -1 1;1 0 0;-1 0 0]; %skew symmetric matrix representing w cross Rwb
z = pi/2
R = expm(Wx*z)