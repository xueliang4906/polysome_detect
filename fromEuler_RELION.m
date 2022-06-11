%% fromEuler TO ROTATION MATRIX

% This is to generate a rotation matrix with 3 euler angles in Relion data star to calculate the converted new xyz values of the input xyz after rotation   

function [A] = fromEuler(rot, tilt, psi)

ca = cosd(rot);
cb = cosd(tilt);
cg = cosd(psi);
sa = sind(rot);
sb = sind(tilt);
sg = sind(psi);
cc = cb * ca;
cs = cb * sa;
sc = sb * ca;
ss = sb * sa;

A = zeros(3);
A(1,1) = cg * cc - sg * sa;
A(1,2) = cg * cs + sg * ca;
A(1,3) = -cg * sb;
A(2,1) = -sg * cc - cg * sa;
A(2,2) = -sg * cs + cg * ca;
A(2,3) = sg * sb;
A(3,1) = sc;
A(3,2) = ss;
A(3,3) = cb;

end

