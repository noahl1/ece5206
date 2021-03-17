function Rx=xrot(phi)

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

%rotation of the magnetization vector about x axis