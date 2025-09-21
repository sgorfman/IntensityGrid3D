function [hgrid,kgrid,lgrid,IG,sigmaIG,VG,CG] = IntensityGrid_plusSymmetry_parallel_Lorentz_V3(FolderName, UB, P0, hL, Nh, kL, Nk, lL, Nl, SymMatrix)
%INTENSITYGRID_PLUSSYMMETRY_PARALLEL_LORENTZ_V3
%   Reconstructs a 3D intensity grid in reciprocal space using symmetry 
%   operations, detector corrections, and Lorentz volume element.
%
% Inputs:
%   FolderName : directory containing .mat image files
%   UB         : 3x3 UB matrix (orientation)
%   P0         : 1xN or 2xN geometry parameter vector
%   hL, kL, lL : [min, max] range for h, k, l
%   Nh, Nk, Nl : number of bins along h, k, l
%   SymMatrix  : 3x3xNS array of symmetry matrices
%
% Outputs:
%   hgrid, kgrid, lgrid : voxel grid edges
%   IG        : average intensity grid
%   sigmaIG   : voxel-wise intensity uncertainty
%   VG        : total voxel volume
%   CG        : number of contributing pixels per voxel

%% 0. Constants and Voxel Grid Initialization

maskPattern = 'UB';                  % File prefix for input .mat files
%mode = 'ID28';                     % Polarization correction mode
mode = 'HomeLab';
detectorSize = [1043, 981];         % PILATUS image shape (rows, cols)

[hgrid, kgrid, lgrid, dh, dk, dl, Nh, Nk, Nl] = prepare_voxel_grid(hL, Nh, kL, Nk, lL, Nl);

[HG, ~, ~] = meshgrid(hgrid, kgrid, lgrid); 
IG = HG .* 0;  % Initialize intensity grid

% Flattened versions
[IG_lin, SG2_lin, Volume_lin, Counter_lin] = deal(zeros(numel(IG),1));
Nbins = numel(IG_lin);

%% 1. Load Hot Pixel Map and Symmetry Info

switch mode
    case 'ID28'
        load('HotPixels_ID28.mat', 'HotPixels_ID28');
        HotPixels = HotPixels_ID28;
    case 'HomeLab'
        load('HotPixels.mat', 'HotPixels');
end

NS = size(SymMatrix,3);

%% 2. Locate All Relevant Data Files

fileList = dir(fullfile(FolderName, '**', [maskPattern, '*.mat']));
Nfiles = numel(fileList);

[X, Y] = meshgrid(1:detectorSize(2), 1:detectorSize(1)); % 981 × 1043

%P = P0(1,:);

hL1 = hL(1) - dh/2; 
kL1 = kL(1) - dk/2; 
lL1 = lL(1) - dl/2; 

%% 3. Begin Parallel Processing Loop

parfor n = 1:Nfiles

    FileName = fullfile(fileList(n).folder, fileList(n).name);
    fprintf('Reading file: %s\n', FileName);

    [I, ~, TTH, Omega, CHI, Phi, DD, DH, dOmega, dPhi, ET] = LoadMATImage(FileName, HotPixels);
    sigmaI2 = abs(I);

    [IG_lin_local, SG2_lin_local, Volume_lin_local, Counter_lin_local] = deal(zeros(Nbins,1,'single'));
    SymMatrix_local = SymMatrix;

    Pn = adjust_geometry(P0,DD,DH,TTH);

    [QX, QY, QZ, TTH_XY, ~, ~, ~, Oblique, ~, ~, Rz] = Pixel_to_XYZ_TAU(X, Y, Pn);

    % Detector correction
    SolidAngle = cosd(Oblique).^3 / Pn(1)^2;
    I = I / ET ./ SolidAngle;
    sigmaI2 = sigmaI2 / ET^2 ./ (SolidAngle.^2);

    % Polarization correction
    switch mode
        case 'ID28'
            PF = (1 - Rz.^2);
        case 'HomeLab'
            PF = (1 + cosd(TTH_XY).^2) / 2;
    end

    I = I ./ PF;
    sigmaI2 = sigmaI2 ./ (PF.^2);

    % Lorentz and geometry matrix
    [Vp, ~, M] = Lorentz(X, Y, Pn, Omega, CHI, Phi, dPhi, dOmega);
    I = I .* Vp;
    sigmaI2 = sigmaI2 .* (Vp.^2);

    % Transform Q to HKL0
    QLAB = [QX(:)'; QY(:)'; QZ(:)'];
    HKL0 = (M * UB) \ QLAB;

    for ns = 1:NS
        Sn = SymMatrix_local(:,:,ns);
        HKL = Sn * HKL0;
        [In, S2n, Vpn] = deal(I, sigmaI2, Vp);

        HN = ceil((HKL(1,:)' - hL1)/dh);
        KN = ceil((HKL(2,:)' - kL1)/dk);
        LN = ceil((HKL(3,:)' - lL1)/dl);

        indexRem = (HN<=0)|(HN>Nh)|(KN<=0)|(KN>Nk)|(LN<=0)|(LN>Nl)|(In(:)<0);
        [HN, KN, LN, In, S2n, Vpn] = deal(HN(~indexRem), KN(~indexRem), LN(~indexRem), In(~indexRem), S2n(~indexRem), Vpn(~indexRem));

        BN = (LN-1)*Nh*Nk + (HN - 1)*Nk + KN;

        for m = 1:numel(BN)
            b = BN(m);
            IG_lin_local(b)  = IG_lin_local(b) + In(m);
            SG2_lin_local(b) = SG2_lin_local(b) + S2n(m);
            Volume_lin_local(b) = Volume_lin_local(b) + Vpn(m);
            Counter_lin_local(b) = Counter_lin_local(b) + 1;
        end    
    end

    IG_lin = IG_lin + IG_lin_local; 
    SG2_lin = SG2_lin + SG2_lin_local;
    Volume_lin = Volume_lin + Volume_lin_local;    
    Counter_lin = Counter_lin + Counter_lin_local;
end

%% 4. Normalize and Reshape Output Grids

indexFilled = (Volume_lin > 0);
IG_lin(indexFilled) = IG_lin(indexFilled) ./ Volume_lin(indexFilled);
SG2_lin(indexFilled) = SG2_lin(indexFilled) ./ (Volume_lin(indexFilled).^2);
SG_lin = sqrt(SG2_lin);

IG = reshape(IG_lin, Nk, Nh, Nl);
sigmaIG = reshape(SG_lin, Nk, Nh, Nl);
VG = reshape(Volume_lin, Nk, Nh, Nl);
CG = reshape(Counter_lin, Nk, Nh, Nl);
