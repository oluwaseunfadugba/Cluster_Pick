function tetra_check
% tetra_check - check the volume calculation for a tetrahedron using a cube

% tetrahedron 1
x1=[0 1 0 0];
y1=[0 0 0 1];
z1=[0 0 1 0];

N=4;
[xn,yn,zn,vol1]=tetvol(x1,y1,z1,N);
fprintf('tetrahedron 1 \n');
fprintf('volume = %g\n',vol1);

% tetrahedron 2
x2=[1 1 0 1];
y2=[0 1 0 0];
z2=[1 1 1 0];

N=4;
[xn,yn,zn,vol2]=tetvol(x2,y2,z2,N);
fprintf('tetrahedron 2 \n');
fprintf('volume = %g\n',vol2);

% tetrahedron 3
x3=[1 1 0 1];
y3=[1 1 1 0];
z3=[0 1 0 0];

N=4;
[xn,yn,zn,vol3]=tetvol(x3,y3,z3,N);
fprintf('tetrahedron 3 \n');
fprintf('volume = %g\n',vol3);

% tetrahedron 4
x4=[0 1 0 0];
y4=[1 1 0 1];
z4=[1 1 1 0];

N=4;
[xn,yn,zn,vol4]=tetvol(x4,y4,z4,N);
fprintf('tetrahedron 4 \n');
fprintf('volume = %g\n',vol4);

% tetrahedron 5 - the inner pyramid
x5=[0 1 0 1];
y5=[1 0 0 1];
z5=[0 0 1 1];

N=4;
[xn,yn,zn,vol5]=tetvol(x5,y5,z5,N);
fprintf('tetrahedron 5 \n');
fprintf('volume = %g\n',vol5);

vol=vol1+vol2+vol3+vol4+vol5;
fprintf('total volume = %g\n',vol);

end

