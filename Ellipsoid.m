function [ nu, rbdy, rint, rext ] = Ellipsoid( A, B, C, N )
%% Ellipsoid.m
%
% Computes the boundary/interior/exterior points and unit outward normals 
% for an ellipsoid using a simple, random distribution method.
%
% Written by A. D. Kim on 4/17/2019

% compute the points on the unit sphere

Y    = randn( N, 3 ) * diag( [ A B C ] );
dist = sqrt( Y(:,1).^2 / A^2 + Y(:,2).^2 / B^2 + Y(:,3).^2 / C^2 );
Y    = Y ./ dist;

% compute the partial derivatives needed for the unit normals

theta   = acos( Y(:,3) / C );
phi     = atan2( Y(:,2) / B, Y(:,1) / A );

x_theta =  A * cos( theta ) .* cos( phi );
y_theta =  B * cos( theta ) .* sin( phi );
z_theta = -C * sin( theta );

x_phi   = -A * sin( theta ) .* sin( phi );
y_phi   =  B * sin( theta ) .* cos( phi );
z_phi   =  0 * cos( theta );
    
% compute the unit normal vectors

Jvec = cross( [ x_theta y_theta z_theta ], [ x_phi y_phi z_phi ], 2 ); 
Jlen = sqrt( Jvec(:,1).^2 + Jvec(:,2).^2 + Jvec(:,3).^2 );
nu   = [ Jvec(:,1)./Jlen, Jvec(:,2)./Jlen, Jvec(:,3)./Jlen ];

% compute interior/exterior points

ell  = 0.25 * min( [ A B C ] );
rbdy = Y;
rint = Y + ell * nu;
rext = Y - ell * nu;

return