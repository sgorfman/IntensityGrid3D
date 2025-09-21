function [hgrid, kgrid, lgrid, dh, dk, dl, Nh, Nk, Nl] = prepare_voxel_grid(hL, Nh, kL, Nk, lL, Nl)
%PREPARE_VOXEL_GRID Creates 1D grid vectors and spacing for reciprocal space voxels.
%
% Inputs:
%   hL - [h_min, h_max] limits for h
%   Nh - number of h divisions
%   kL - [k_min, k_max] limits for k
%   Nk - number of k divisions
%   lL - [l_min, l_max] limits for l
%   Nl - number of l divisions
%
% Outputs:
%   hgrid - grid vector for h (length Nh+1)
%   kgrid - grid vector for k (length Nk+1)
%   lgrid - grid vector for l (length Nl+1)
%   dh, dk, dl - spacing between voxels along h, k, l
%   Nh, Nk, Nl - incremented grid counts (for edge inclusion)

    % Create linearly spaced grid edges
    hgrid = linspace(hL(1), hL(2), Nh + 1);
    kgrid = linspace(kL(1), kL(2), Nk + 1);
    lgrid = linspace(lL(1), lL(2), Nl + 1);

    % Calculate voxel size (uniform spacing assumed)
    dh = hgrid(2) - hgrid(1);
    dk = kgrid(2) - kgrid(1);
    dl = lgrid(2) - lgrid(1);

    % Update counts to include all edges
    Nh = Nh + 1;
    Nk = Nk + 1;
    Nl = Nl + 1;
end
