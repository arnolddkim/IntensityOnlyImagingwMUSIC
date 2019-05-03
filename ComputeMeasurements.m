function Udata = ComputeMeasurements( k0, csca, rext, X_Rx, Y_Rx, Z_Rx, ...
    Nobjects )
%% ComputeMeasurements.m
%
% Compute the scattered field on the receiver plane by evaluating the 
% Method of Fundamental Solutions approximation.
%
% Written by A. D. Kim on 4/17/2019

% compute useful index arrays

[ indx, jndx ] = meshgrid( ( 1 : length( X_Rx ) ), ( 1 : length( csca ) ) );

% allocate memory for the measurements

Udata = zeros( length( X_Rx ), 1 );

% loop over all objects

for nobjs = 1 : Nobjects

    % extract MFS exterior point data into coordinates

    xext = rext(:,1,nobjs);
    yext = rext(:,2,nobjs);
    zext = rext(:,3,nobjs);

    % compute the distance from Rx to MFS exterior points

    R = sqrt( ( X_Rx(indx) - xext(jndx) ).^2 ...
        + ( Y_Rx(indx) - yext(jndx) ).^2 ...
        + ( Z_Rx - zext(jndx) ).^2 );

    % compute Green's function for these distances

    G = exp( 1i * k0 * R ) ./ ( 4 * pi * R );

    % compute the scattered field on the Rx

    Udata = Udata + G.' * csca(:,nobjs);
    
end

return