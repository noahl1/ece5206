function Rth=throt(phi,theta)
%returns the rotation matrix for a rotation of phi about any axis
Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;