function [Data, Header, TTH, Omega, Chi, Phi, det_d, det_vo, dOmega, dPhi, ET] = LoadMATImage(ImageName, HotPixels, Frame)
%LOADMATIMAGE Loads detector image and metadata, applies hot pixel correction,
% and extracts relevant geometric parameters including scan increments.
%
% Inputs:
%   ImageName  - full path to .mat file
%   HotPixels  - mask or list for hot pixel removal
%   Frame      - (optional) 2x2 ROI window [Xstart Xend; Ystart Yend]
%
% Outputs:
%   Data       - cleaned image data
%   Header     - detector metadata
%   TTH        - 2-theta angle
%   Omega      - sample rotation angle
%   Chi        - sample tilt
%   Phi        - azimuthal angle
%   det_d      - detector distance
%   det_vo     - detector vertical offset
%   dPhi       - Phi increment (rotation step)
%   dOmega     - Omega increment (rotation step)

    %% Load image and header
    loaded = load(ImageName, 'F');
    Data = loaded.F.data;
    Header = loaded.F.header;

    %% Reorient and clean data
    Data = fliplr(Data)';                    % Rotate image for display
    Data = RemoveHotPixels(Data, HotPixels); % Apply hot pixel correction

    %% Optional cropping
    if nargin > 2 && ~isempty(Frame)
        Data = Data(Frame(2,1):Frame(2,2), Frame(1,1):Frame(1,2));
    end

    %% Metadata extraction with safe defaults
    TTH     = safe_hdr_extract(Header, 'Detector_2theta ', 0);
    Omega   = safe_hdr_extract(Header, 'Omega ', 0);
    Chi     = safe_hdr_extract(Header, 'Chi ', 0);
    Phi     = safe_hdr_extract(Header, 'Phi ', 0);
    det_d   = safe_hdr_extract(Header, 'Detector_distance ', NaN);
    det_vo  = safe_hdr_extract(Header, 'Detector_Voffset ', NaN);
    dPhi    = safe_hdr_extract(Header, 'Phi_increment ', 0);      % default increment
    dOmega  = safe_hdr_extract(Header, 'Omega_increment ', 0);    % default increment
    ET      = safe_hdr_extract(Header, 'Exposure_time ', 0);      % default increment

end
