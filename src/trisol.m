function [alph,beta,gama]=trisol(a,b,c,flag)

% Solve for the angles of a triangle given the length of each side

s=0.5.*(a+b+c);

% To avoid near zero negative values. 
% Happens when an earthquake is on an edge of a fault
sa = real(abs(s-a));
sb = real(abs(s-b));
sc = real(abs(s-c));

r=real(sqrt(sa.*sb.*sc./s));

alph=2.0.*atan2(r,sa);
beta=2.0.*atan2(r,sb);
gama=2.0.*atan2(r,sc);

if flag == 'd'
    con=180./pi;
    alph=alph.*con;
    beta=beta.*con;
    gama=gama.*con;
end
return;