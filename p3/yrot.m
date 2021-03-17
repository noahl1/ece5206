function Ry=yrot(phi)

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];

%rotation of the magnetization vector about y axis