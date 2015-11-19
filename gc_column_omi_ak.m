function [ gc_no2 ] = gc_column_omi_ak( gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp )
%GC_COLUMN_OMI_AK Applies OMI averaging kernels to calculation of GEOS-Chem columns.
%   Averaging kernels are used to describe how a change in the "true" state
%   of a system is expressed as a change in the observations of a system.
%   In terms of satellite observations, a modeled NO2 column is usually
%   taken as "truth" (whether it is actually true is another matter
%   entirely) and the satellite column as observations. Thus, a vector of
%   averaging kernels is used to weight a modeled column for an
%   apples-to-apples comparison with a satellite observation.
%
%   For satellite retrievals, averaging kernels are defined as a vector of
%   box air mass factors or scattering weights scaled by the inverse of the
%   total AMF.  From Burrows, Platt, and Borrel (ed.), this is specifically
%   defined as:
%
%       AK_i = BAMF_i / AMF
%
%   and the box air mass factor (BAMF) is at one point defined as:
%
%       BAMF_i = SCD_i / VCD_i
%
%   therefore is the air mass factor for level i only. The observation VCD
%   can then be calculated from the true VCD as:
%
%       VCD_obs = sum_i(AK_i * VCD_i) = 1/AMF * sum_i(BAMF_i * VCD_i)
%
%   This is essentially calculating the observed slant column (the result
%   of BAMF_i * VCD_i summed over all levels) then dividing by the AMF to
%   simulate the retrieval process.
%
%   For OMI, the box air mass factors are given instead as scattering
%   weights, which reflect a similar quantity: basically the sensitivity of
%   the satellite to NO2 at that pressure level. From Eskes and Boersma
%   (2003), we know that the BAMF is independent of trace gas distribution
%   for optically thin absorbers, and NO2 is an optically thin absorber in
%   the atmosphere. Therefore, the scattering weights (which depend on
%   optical properties of the atmosphere & earth plus the geometric
%   arrangement of the sun and satellite) are sufficient for this purpose.
%
%   Reading:
%       Burrows, John P.; Platt, Ulrich; Borrell, Peter (eds.). The Remote
%       Sensing of Tropospheric Composition from Space. Springer, 2011. see
%       pp. 91-97.
%
%       Eskes, H.J. and Boersma, K.F. Averaging Kernels for DOAS
%       total-column satellite retrievals. Atmos. Chem. Phys. (2003), 3,
%       1285-1291.
%
%
%   This function will apply OMNO2 averaging kernels to GEOS-Chem output.
%   It requires as input 4 output structures from GEOS-Chem created using
%   read_geos_output or by using read_gc_nd51.py and
%   convert_py_nd51_structure.m: the NO2 mixing ratios, the matrix of
%   pressures, the number density of air, and the tropopause level. It will
%   assume that you want to apply the averaging kernels to the full time
%   extent of the NO2 structure.
%
%   For each day, this function will load all 14 OMI orbits then average
%   the vectors of averaging kernels to the GEOS-Chem grid boxes.  The
%   averaging kernels 
%
%   Dependencies:
%       Classes/JLLErrors.m

% Get run directories
global onCluster
if isempty(onCluster)
    onCluster = false;
end

fprintf('onCluster = %d\n',onCluster);

global omi_he5_dir
if onCluster && isempty(omi_he5_dir)
    E.runscript_error('omi_he5_dir')
elseif ~onCluster
    omi_he5_dir = '/Volumes/share-sat/SAT/OMI/OMNO2';
end

fprintf('OMI dir = %s\n',omi_he5_dir');

[gc_loncorn, gc_latcorn] = geos_chem_corners;

tVec = gc_no2.tVec;
gc_no2_data = gc_no2.dataBlock;
gc_bxhght_data = gc_bxhght.dataBlock;
gc_pressure_data = gc_pressure.dataBlock;
gc_ndens_air_data = gc_ndens_air.dataBlock;
gc_tp_data = gc_tp.dataBlock;
columns = nan(size(gc_tp_data));

parfor d=1:numel(tVec)
    fprintf('Loading OMI files for %s\n',datestr(tVec(d)));
    [omi_aks, omi_lon, omi_lat] = load_omi_files(year(tVec(d)), month(tVec(d)), day(tVec(d)), omi_he5_dir);
    fprintf('Binnind OMI AKs for %s\n',datestr(tVec(d)));
    binned_aks = bin_omi_aks(gc_loncorn, gc_latcorn, omi_aks, omi_lon, omi_lat);

    % binned_aks is output in ak_vec x gc_nlat x gc_nlon x time, rearrange to
    % gc_nlon x gc_nlat x ak_vec x time.
    binned_aks = permute(binned_aks, [3 2 1 4]);
    
    %tmp = integrate_geoschem_profile(cat(1,gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp), 1e-9, [], 0, binned_aks);
    %gc_no2.Columns = tmp(1).Columns;
    
    columns(:,:,d) = integrate_omi_profile(gc_no2_data(:,:,:,d), gc_bxhght_data(:,:,:,d), gc_pressure_data(:,:,:,d), gc_ndens_air_data(:,:,:,d), gc_tp_data(:,:,d), binned_aks);
end

gc_no2.Columns = columns;

end

function [omi_aks, omi_lon, omi_lat] = load_omi_files(yr, mn, dy, omi_he5_dir)
E = JLLErrors;

full_path = fullfile(omi_he5_dir,sprintf('%04d',yr),sprintf('%02d',mn));
file_pattern = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',yr,mn,dy);
F = dir(fullfile(full_path, file_pattern));

fprintf('%s\n',full_path);
fprintf('%s\n',file_pattern);

if isempty(F)
    E.filenotfound(sprintf('OMI files for %04d-%02d-%02d',yr,mn,dy));
end


% We'll concatenate along the along-track dimension since that varys with
% each swath and the other two are constant.
omi_aks = [];
omi_lon = [];
omi_lat = [];

for a=1:numel(F)
    %fprintf('Loading file %d of %d\n',a,numel(F));
    hi = h5info(fullfile(full_path, F(a).name));
    omi.sw = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'ScatteringWeight')));
    omi.sw(omi.sw < 1e-29) = nan;
    omi.amf = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AmfTrop')));
    omi.amf(omi.amf < 1e-29) = nan;
    
    % These are fields needed to reject pixels
    omi.CloudFraction = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'CloudFraction')))*1e-3;
    omi.CloudFraction(omi.CloudFraction<0) = nan;
    omi.ColumnAmountNO2Trop = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'ColumnAmountNO2Trop')));
    omi.ColumnAmountNO2Trop(omi.ColumnAmountNO2Trop<-1e29) = nan;
    omi.vcdQualityFlags = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1, 'VcdQualityFlags')));
    omi.vcdQualityFlags(omi.vcdQualityFlags==65535)=nan;
    omi.XTrackQualityFlags = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'XTrackQualityFlags')));
    omi.XTrackQualityFlags(omi.XTrackQualityFlags==255) = nan;
    omi.Areaweight = ones(size(omi.amf));
    
    omi = omi_sp_pixel_reject(omi,'geo',0.2,'XTrackFlags');
    
    rejects = omi.Areaweight == 0;
    omi.sw(:,rejects) = nan;
    omi.amf(rejects) = nan;
    
    omi_sw = omi.sw;
    omi_amf = omi.amf;
    
    % Remove the large 'omi' structure since we have gotten what we need
    % out of it.
    clear('omi');
    
    this_aks = nan(size(omi_sw));
    for x=1:size(omi_amf,1)
        for y=1:size(omi_amf,2)
            this_aks(:,x,y) = omi_sw(:,x,y) ./ omi_amf(x,y);
        end
    end
    
    omi_aks = cat(3, omi_aks, this_aks);
    
    this_lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
    this_lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
    omi_lon = cat(2, omi_lon, this_lon);
    omi_lat = cat(2, omi_lat, this_lat);
end

end

function aks = bin_omi_aks(gc_loncorn, gc_latcorn, omi_aks, omi_lon, omi_lat)
% This subfunction will handle the binning of OMNO2 averaging kernels to
% the GEOS-Chem grid. Inputs: matrices of GEOS-Chem corner points and
% matrices of OMI pixel AKs (with the different AK levels along the first
% dimension) and pixel lon/lats.  gc_loncorn and gc_latcorn are expected to
% be matrices that are (m+1)-by-(n+1) if there are m-by-n GC grid cells.

E = JLLErrors;

if any(size(gc_loncorn) ~= size(gc_latcorn))
    E.sizeMismatch('gc_loncorn','gc_latcorn');
end

sz_aks = size(omi_aks);
sz_omilon = size(omi_lon);
sz_omilat = size(omi_lat);

if ndims(omi_aks) ~= 3
    E.badinput('omi_aks expected to be 3 dimensional')
elseif ~ismatrix(omi_lon)
    E.badinput('omi_lon expected to be 2 dimensional')
elseif ~ismatrix(omi_lat)
    E.badinput('omi_lat expected to be 2 dimensional')
end

% There should usually be fewer points in the latitude direction than the
% longitude direction of GEOS-Chem, use this to make sure the matrices have
% the right orientation.
if size(gc_loncorn,1) > size(gc_loncorn,2)
    E.badinput('The GC corner matrices are expected to be lat x lon')
end
if any(diff(gc_loncorn(1,:)) < 0) || any(diff(gc_latcorn(:,1)) < 0)
    E.badinput('The GC corner matrices are expected to have the SW corner at (1,1) and NE corner at (end,end)')
end


if any(sz_aks(2:3)~=sz_omilon) || any(sz_aks(2:3)~=sz_omilat)
    E.badinput('The second and third dimensions of omi_aks should be the same as the dimensions of omi_lon and omi_lat')
end

% Loop over all GEOS-Chem cells and bin the averaging kernels to them.
gc_nlat = size(gc_loncorn,1)-1;
gc_nlon = size(gc_loncorn,2)-1;
aks = nan(size(omi_aks,1), gc_nlat, gc_nlon);
for a=1:gc_nlat
    %fprintf('Binning %.1f%% complete\n',a/gc_nlat*100);
    for b=1:gc_nlon
        x1 = gc_loncorn(a,b);
        x2 = gc_loncorn(a,b+1);
        y1 = gc_latcorn(a,b);
        y2 = gc_latcorn(a+1,b);
        
        xx = omi_lon >= x1 & omi_lon < x2 & omi_lat >= y1 & omi_lat < y2;
        % All OMI aks are given at the same pressures, so we can just
        % average across all kernels that belong to this grid cell.
        this_ak = nanmean(omi_aks(:,xx),2);
        aks(:,a,b) = this_ak;
    end
end


end

function no2_columns = integrate_omi_profile(gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp, binned_aks)
% This will integrate the geos-chem columns weighted by the OMI AKs after
% interpolating to the OMI pressures.  The integration will be carried out
% after that in the same manner as integrate_geoschem_profile, that is,
% the integral will be calculated using the rectangular rule.
%
% In order to find the box height in meters, the height of each GC pressure
% level will be calculated as the cumulative sum of the box height values,
% which can then be intepolated to the OMI pressure levels, and the box
% height for level i taken as the difference between level i and i+1
% altitude in meters. Number density will also be interpolated to the new
% pressure levels, and the GC tropopause will be used to restrict the
% column.

% First, handle all the necessary intepolation to regrid the GC variables
% to the OMI pressure levels.
omi_pres = [1020 1010 1000 990 975 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200 120 60 35 20 12 8 5 3 1.5 0.8];
omi_pres_edge = [1025 1015 1005 995 982.5 967.5 952.5 935 912.5 887.5 862.5 837.5 812.5 785 755 720 680 635 585 525 475 425 375 315 240 150 72.5 42.5 24 14 9.5 6 3.75 1.85 1.15 0.45];

sz = size(gc_no2);
if numel(sz) < 4;
    sz = cat(2, sz, ones(1,4-numel(sz)));
end
gc_z = cumsum(gc_bxhght,3);
omi_bxhght = nan(sz(1), sz(2), numel(omi_pres), sz(4));
omi_ndens_air = nan(sz(1), sz(2), numel(omi_pres), sz(4));
omi_no2 = nan(sz(1), sz(2), numel(omi_pres), sz(4));
for a=1:sz(1)
    for b=1:sz(2)
        for t=1:sz(4)
            % For each column, interpolate box altitude to OMI pressures,
            % height is proportional to ln(p), so use linear interpolation
            % w.r.t. ln(p):
            %   z = -ln(p/p0)*H = H*log(p0) - H*log(p)
            % When GC outputs actual pressure edges the last edge is 0.01
            % which is not included in the satellite pressure values, we'll
            % include it to have edges to work with, then also add 0 to the
            % boxheight values so that we can take the difference.
            gc_pvec = cat(1,squeeze(gc_pressure(a,b,:,t)),0.01);
            gc_zvec = cat(1, 0, squeeze(gc_z(a,b,:,t)));
            omi_zvec = interp1(log(gc_pvec), gc_zvec, log(omi_pres_edge),'linear','extrap');
            omi_bxhght(a,b,:,t) = diff(omi_zvec);
            
            % Also interpolate number density to OMI pressure levels,
            % number density is directly proportional to p:
            %   pV = nRT => n/V = p/RT
            % so we'll directly interpolate.
            gc_nair_vec = squeeze(gc_ndens_air(a,b,:,t));
            omi_ndens_air(a,b,:,t) = interp1(gc_pvec(1:end-1), gc_nair_vec, omi_pres, 'linear','extrap');
            
            % Finally interpolate the GC NO2 mixing ratio values (scaling
            % to straight mixing ratios, i.e. parts-per-part). I'll do this
            % in log-log space, which I've seen before. We will NOT
            % extrapolate this time, because I do not want there to be
            % values below the surface (so outside of the GC pressure
            % values, it will have the value of NaN)
            gc_no2vec = squeeze(gc_no2(a,b,:,t));
            omi_no2(a,b,:,t) = exp(interp1(log(gc_pvec(1:end-1)), log(gc_no2vec), log(omi_pres)));
            strat = omi_pres < gc_pressure(a,b,floor(gc_tp(a,b,t)),t);
            omi_no2(a,b,strat,t) = nan;
        end
    end
end

% Finally the actual calculation of columns is fairly straightforward: sum
% over each box's partial column weighted by AK (after appropriate unit
% conversions)

omi_no2_ndens = (omi_no2 * 1e-9) .* (omi_ndens_air * 1e-6) .* (omi_bxhght * 100) .* binned_aks;
no2_columns = squeeze(nansum2(omi_no2_ndens,3));

end
