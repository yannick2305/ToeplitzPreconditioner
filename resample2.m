function pointsUni = resample2(points,m)
% Resample complex points by uniform angle

%points

x = real(points);
y = imag(points);
k=3;
r = sampleUniformlyVec(x,y,m,k);

thetaUni = linspace(0,pi,m);

r2 = [r([m-1:-1:2]),r];
theta2 = [2*pi-thetaUni([m-1:-1:2]),thetaUni];

xUni = r2 .* cos(theta2);
yUni = r2 .* sin(theta2);

pointsUni = complex(xUni,yUni);


