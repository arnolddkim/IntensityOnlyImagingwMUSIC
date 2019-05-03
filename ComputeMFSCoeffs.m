function [ cint, csca ] = ComputeMFSCoeffs( rbdy, rint, rext, nu, dir, ...
    k0, k1, N, Nobjects )
%% ComputeMFSCoeffs.m
%
% Compute the solution of the linear system of equations that give the MFS
% expansion coefficients resulting from a collocation method to approximate
% the boundary conditions.
%
% Written by A. D. Kim on 4/17/2019

% set the size of the linear system to be solved

NN = 2 * N * Nobjects;

% allocate memory for the linear system to be solved

bvec = zeros( NN, 1 );
Amtx = zeros( NN, NN );

% compute useful index arrays

[ indx, jndx ] = meshgrid( ( 1 : N ) );

% loop over all scattering objects

for mobjs = 1 : Nobjects

    % compute index arrays
    
    indx1 = 2 * N * ( mobjs - 1 ) + ( 1 : N );
    indx2 = 2 * N * ( mobjs - 1 ) + ( N + 1 : 2 * N );

    % compute the incident field propagating in direction dir

    dirDOTrbdy = dir(:,1) .* rbdy(:,1,mobjs) ...
        + dir(:,2) .* rbdy(:,2,mobjs) ...
        + dir(:,3) .* rbdy(:,3,mobjs);

    dirDOTnu = dir(:,1) .* nu(:,1,mobjs) ...
        + dir(:,2) .* nu(:,2,mobjs) ...
        + dir(:,3) .* nu(:,3,mobjs);

    bvec(indx1) = exp( 1i * k0 * dirDOTrbdy );
    bvec(indx2) = 1i * k0 * dirDOTnu .* exp( 1i * k0 * dirDOTrbdy );

    % --------------------------------- %
    % blocks of Amtx for interior field %
    % --------------------------------- %
    
    % compute the interior distances and cosines

    Rint = reshape( sqrt( ( rbdy(indx,1,mobjs) - rint(jndx,1,mobjs) ).^2 ...
        + ( rbdy(indx,2,mobjs) - rint(jndx,2,mobjs) ).^2 ...
        + ( rbdy(indx,3,mobjs) - rint(jndx,3,mobjs) ).^2 ), N, N );

    MUint = reshape( nu(indx,1,mobjs) ...
        .* ( rbdy(indx,1,mobjs) - rint(jndx,1,mobjs) ) ...
        + nu(indx,2,mobjs) ...
        .* ( rbdy(indx,2,mobjs) - rint(jndx,2,mobjs) ) ...
        + nu(indx,3,mobjs) ...
        .* ( rbdy(indx,3,mobjs) - rint(jndx,3,mobjs) ), N, N ) ./ Rint;

    % assign the entries of the block

    Amtx(indx1,indx2) = ComputeG( k1, Rint );
    Amtx(indx2,indx2) = ComputeDnG( k1, Rint, MUint );           

    % --------------------------------- %
    % blocks of Amtx for exterior field %
    % --------------------------------- %

    for nobjs = 1 : Nobjects
    
        % compute index arrays
        
        indx3 = 2 * N * ( nobjs - 1 ) + ( 1 : N );

        % compute the exterior distances and cosines

        Rext = reshape( sqrt( ( rbdy(indx,1,mobjs) - rext(jndx,1,nobjs) ).^2 ...
            + ( rbdy(indx,2,mobjs) - rext(jndx,2,nobjs) ).^2 ...
            + ( rbdy(indx,3,mobjs) - rext(jndx,3,nobjs) ).^2 ), N, N );

        MUext = reshape( nu(indx,1,mobjs) ...
            .* ( rbdy(indx,1,mobjs) - rext(jndx,1,nobjs) ) ...
            + nu(indx,2,mobjs) ...
            .* ( rbdy(indx,2,mobjs) - rext(jndx,2,nobjs) ) ...
            + nu(indx,3,mobjs) ...
            .* ( rbdy(indx,3,mobjs) - rext(jndx,3,nobjs) ), N, N ) ./ Rext;
        
        % assign the entries of the block matrices
        
        Amtx(indx1,indx3) = -ComputeG( k0, Rext );
        Amtx(indx2,indx3) = -ComputeDnG( k0, Rext, MUext );
        
    end
    
end

%% solve the linear system for the expansion coefficients

c = Amtx \ bvec;

%% re-organize the expansion coefficients

csca = zeros( N, Nobjects );
cint = zeros( N, Nobjects );

for mobjs = 1 : Nobjects
    
    indx1 = 2 * N * ( mobjs - 1 ) + ( 1 : N );
    indx2 = 2 * N * ( mobjs - 1 ) + ( N + 1 : 2 * N );
   
    csca(:,mobjs) = c(indx1);
    cint(:,mobjs) = c(indx2);
    
end

return

function G = ComputeG( k, R )
%% ComputeG.m
%
% Compute the whole space Green's function for wavenumber k at radial
% distance R.
%
% Written by A. D. Kim on 4/17/2019

G = exp( 1i * k * R ) ./ ( 4 * pi * R );

return

function DnG = ComputeDnG( k, R, CosTheta )
%% ComputeG.m
%
% Compute the whole space Green's function for wavenumber k at radial
% distance R.
%
% Written by A. D. Kim on 4/17/2019

% compute Green's function

G = ComputeG( k, R );

% compute the normal derivative of Green's function

DnG = CosTheta .* ( 1i * k - 1 ./ R ) .* G;

return