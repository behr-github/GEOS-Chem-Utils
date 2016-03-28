function [ varargout ] = misc_gc_plotting_fxns( plttype, varargin )
%GC_PLOTTING_FXNS Misc. plots to make for GC outputs
%   Rather than trying to do this all in the command window, or having a
%   zillion different .m files all over the place, I'm just going to
%   collect them all here.
%
%   Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2015

DEBUG_LEVEL = 2;
E = JLLErrors;

plttype = lower(plttype);
switch plttype
    case 'subfrac_lnox'
        [varargout{1}, varargout{2}] = subfrac_lnox(varargin{1});
    case 'trop_mask'
        [varargout{1}] = trop_mask(varargin{1:2});
    case 'region_xxyy'
        [varargout{1}, varargout{2}] = region_xxyy(varargin{1});
        [~,~,varargout{3}] = define_regions(varargin{1});
    case 'timeser_column'
        plot_dcolumn_percentiles(varargin{:}, 'relative', false);
    case 'timeser_relcol'
        plot_dcolumn_percentiles(varargin{:}, 'relative', true);
    case 'timeser_plev'
        plot_level_dconc(varargin{:}, 'relative', false);
    case 'timeser_relplev'
        plot_level_dconc(varargin{:}, 'relative', true);
    case 'timeser_dlnox'
        plot_lnox_col_enhnc_diff(varargin{:}, 'relative', false);
    case 'timeser_reldlnox'
        plot_lnox_col_enhnc_diff(varargin{:}, 'relative', true);
    case 'meridional'
        plot_meridional_avg(varargin{:});
    case 'meridional-sat'
        add_sat_mer_avg_to_fig(varargin{:});
    case 'meridional_table'
        varargout{1} = avg_meridional_avg_table(varargin{:});
    case 'grid_mls'
        varargout{1} = grid_mls_data_to_gc(varargin{:});
    case 'grid_mls_daily'
        [varargout{1}, varargout{2}] = grid_daily_mls_data_to_gc(varargin{:});
    case 'grid_month_omno2'
        [varargout{1}, varargout{2}] = grid_omno2d_monthly(varargin{:});
    case 'grid_day_omno2'
        varargout{1} = grid_omno2d_daily(varargin{:});
    case 'comp2sat'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}] = compare_regions_to_sat(varargin{:});
    otherwise
        fprintf('Did not recognize plot type\n');
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OTHER FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [prodfrac, lnoxgt60] = subfrac_lnox(prodfrac)
        % Calculates the fraction of NO emission due to lightning only
        % considering anthropogenic and biomass burning as other sources
        % Takes a production structure output from geos_frac_lnox
        subtot = prodfrac.NO_molec_cm2_sec.Anthropogenic_NO + prodfrac.NO_molec_cm2_sec.Biomass_NO + prodfrac.NO_molec_cm2_sec.Lightning_NO;
        subfrac = prodfrac.NO_molec_cm2_sec.Lightning_NO ./ subtot;
        prodfrac.NO_frac.Subfrac_LNOx = subfrac;
        lnoxgt60 = subfrac > 0.6;
    end

    function dataind = get_data_inds(casestruct)
        % Finds which indices correspond to pressure surface, box height,
        % number density, and tropopause level. Return as a struct with
        % fields bxhght, psurf, ndens, and tplev.
        dataind = struct('bxhght',nan,'psurf',nan,'ndens',nan,'tplev',nan);
        for n=1:numel(casestruct)
            if ~isempty(regexpi(casestruct(n).fullName,'bxheight'));
                dataind.bxhght = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'psurf|pressure'))
                dataind.psurf = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'number density'))
                dataind.ndens = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'tropopause'))
                dataind.tplev = n;
            end
        end
    end

    function [relative, ndens_bool, lonlim, latlim, timelim_bool] = parse_varargs(vargs)
        % Parse varargs. If one is the string 'relative', the next one is
        % the relative boolean. Also look for the string 'ndens', if there,
        % set it to use number density (and remove it). There should be 0
        % to 2 left after that, if 2, then those are the lat and lon
        % limits. If 1, it should be a region definition
        xx = strcmpi('relative',vargs);
        if sum(xx) == 1
            ff = find(xx);
            relative = vargs{ff+1};
            vargs(ff:(ff+1)) = [];
        elseif sum(xx) == 0
            relative = false;
        else
            E.badinput('The relative keyword appears more than once');
        end
        
        xx = strcmpi('ndens',vargs);
        if sum(xx) > 0
            ndens_bool = true;
            vargs(xx) = [];
        elseif sum(xx) == 0
            ndens_bool = false;
        end
        
        xx = strcmpi('timelim',vargs);
        if sum(xx) > 0
            timelim_bool = true;
        else
            timelim_bool = false;
        end
        
        if numel(vargs) == 0
            lonlim = [-200 200];
            latlim = [-100 100];
        elseif numel(vargs) == 2
            lonlim = vargs{1};
            latlim = vargs{2};
        elseif numel(vargs) == 1 && ischar(vargs{1}) && ismember(vargs{1}, {'na','sa','naf','saf','seas','neur'});
            [lonlim, latlim] = define_regions(vargs{1});
        else
            E.badinput('Parsing of optional arguments failed')
        end
    end

    function [lonlim, latlim, timeind] = define_regions(region)
        % If limiting by time, define the Northern Hemisphere to be
        % June-Aug and the Southern Hemisphere to be Jan, Feb, Dec.
        nh_times = modis_date_to_day('2012-06-01'):modis_date_to_day('2012-08-31');
        nh_mask = false(1,366);
        nh_mask(nh_times) = true;
        sh_times = [modis_date_to_day('2012-01-01'):modis_date_to_day('2012-02-29'), modis_date_to_day('2012-12-01'):modis_date_to_day('2012-12-31')];
        sh_mask = false(1,366);
        sh_mask(sh_times) = true;
        switch region
            case 'na'
                lonlim = [-120, -65];
                latlim = [20, 60];
                %latlim = [20, 50];
                timeind = nh_mask;
            case 'sa'
                lonlim = [-77, -39];
                latlim = [-35, 10];
                timeind = sh_mask;
            case 'naf'
                lonlim = [-15, 48];
                latlim = [3, 25];
                timeind = nh_mask;
            case 'saf'
                lonlim = [10, 48];
                latlim = [-30, 3];
                timeind = sh_mask;
            case 'seas'
                lonlim = [95, 146];
                latlim = [-9, 26];
                timeind = nh_mask;
            case 'neur'
                lonlim = [60, 130];
                latlim = [30, 68];
                timeind = nh_mask;
            case 'nh'
                lonlim = [-180,180];
                latlim = [0, 60];
                timeind = nh_mask;
            case 'sh'
                lonlim = [-180, 180];
                latlim = [-60, 60];
                timeind = sh_mask;
            otherwise
                lonlim = [-200, 200];
                latlim = [-100, 100];
                timeind = true(1,366);
        end
    end

    function [xx,yy] = region_xxyy(region)
        [lonlim, latlim] = define_regions(region);
        [glon, glat] = geos_chem_centers('2x25');
        xx = glon >= lonlim(1) & glon <= lonlim(2);
        yy = glat >= latlim(1) & glat <= latlim(2);
    end

    function [mask] = trop_mask(sz, TP)
        % Returns a mask (logical matrix) that is 1 for any grid box in the
        % troposphere (below the tropopause level) and 0 otherwise.
        % Requires the size of the mask desired (should match the size of
        % the dataBlock it will be applied to) and TP is a structure output
        % from read_geos_output with the tropopause level as the dataBlock.
        % The first, second, and final dimensions of the TP dataBlock must
        % be the same size as the corresponding dimensions in the mask
        % desired.
        
        if ~isvector(sz) || numel(sz) < 3 || numel(sz) > 4
            E.badinput('sz must be a size vector with 3 or 4 elements')
        elseif ~isstruct(TP) || ~isfield(TP,'dataBlock') || ~isfield(TP,'fullName') || ~isscalar(TP) 
            E.badinput('TP must be a scalar structure output from read_geos_output')
        elseif isempty(regexpi(TP.fullName,'tropopause'))
            E.badinput('TP does not seem to contain the tropopause level output')
        end
        
        sz_tp = size(TP.dataBlock);
        if numel(sz_tp) > 3 || numel(sz_tp) < 2
            E.badinput('The tropopause dataBlock is expected to have 2 or 3 dimensions');
        elseif any(sz_tp(1:2) ~= sz(1:2)) 
            E.badinput('The first two dimensions of the requested mask size must be the same as the first two dimensions of the TP dataBlock')
        elseif sz_tp(3) ~= sz(end)
            E.badinput('The final dimension of the requested mask is expected to be the same size as the third dimension of the TP dataBlock')
        end
        
        tplev = floor(TP.dataBlock);
        
        mask = true(sz);
        
        for t=1:sz_tp(3)
            for a=1:sz_tp(1)
                for b=1:sz_tp(2)
                    mask(a,b,tplev(a,b,t):end,t) = false;
                end
            end
        end
    end

    

    function [sat_lon, sat_avg] = sat_meridional_avg(specie, region, gridded_data, lnox_crit, noregrid_bool, timelim_bool)
        % Calculates a meridional average for satellite measurements of
        % NO2 (column), O3, or HNO3. Needs which species (as a string),
        % one of the regions (also a string). Hard coded to do 2012 right
        % now.
        [xx,yy] = region_xxyy(region);
        % Add one more point because we actually want to subset the corners
        % not the centers.
        xx_centers = xx; % we'll need this at the end.
        yy_centers = yy;
        xx(find(xx,1,'last')+1) = true;
        yy(find(yy,1,'last')+1) = true;
        [gloncorn, glatcorn] = geos_chem_corners;
        % geos_chem_corners outputs in lat x lon, so transpose just because
        % I'm more used (now) to having lon first.
        gloncorn = gloncorn(yy,xx)';
        glatcorn = glatcorn(yy,xx)';
        
        [lonlim, latlim, timemask] = define_regions(region);
        
        switch lower(specie)
            case 'no2'
                if ~exist('gridded_data','var') || isempty(gridded_data)
                    addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/OMI Utils');
                    gridded_data = omno2d_timeavg('2012-01-01','2012-12-31');
                end
                [omno2d_lon, omno2d_lat] = omno2d_centers;
                
                xx_omi = omno2d_lon(:,1) >= lonlim(1) & omno2d_lon(:,1) <= lonlim(2);
                yy_omi = omno2d_lat(1,:) >= latlim(1) & omno2d_lat(1,:) <= latlim(2);
                
                if exist('lnox_crit','var') && ~isempty(lnox_crit)
                    gridded_data(~lnox_crit) = nan;
                end
                
                if timelim_bool
                    if size(gridded_data,3) == 366
                        gridded_data(:,:,~timemask) = nan;
                    else
                        E.badinput('If you wish to filter by time, the input data must represent 366 days');
                    end
                end
                
                if size(gridded_data,3)>1
                    %gridded_data = nanmean(gridded_data,3);
                    gridded_data = nanmedian(gridded_data,3);
                end
                
                if size(gridded_data,1) * size(gridded_data,2) > 144*91
                    if noregrid_bool
                        gridded_data = gridded_data(xx_omi,yy_omi);
                    else
                        gridded_data = grid_omno2d_to_gc(gridded_data, omno2d_lon, omno2d_lat, gloncorn, glatcorn);
                    end
                else
                    gridded_data = gridded_data(xx_centers, yy_centers);
                end

                if noregrid_bool
                    sat_lon = omno2d_lon(xx_omi,1);
                else
                    sat_lon = geos_chem_centers('2x25');
                    sat_lon = sat_lon(xx_centers);
                end
            otherwise
                if ~exist('gridded_data','var') || isempty(gridded_data)
                    gridded_data = grid_mls_data_to_gc(gloncorn, glatcorn,specie);
                    if timelim_bool
                        warning('Time limits not applied; you must pass an existing matrix for that to work')
                    end
                else
                    if exist('lnox_crit','var') && ~isempty(lnox_crit)
                        gridded_data(~lnox_crit) = nan;
                    end
                    
                    if timelim_bool
                        if size(gridded_data,3) == 366
                            gridded_data(:,:,~timemask) = nan;
                        else
                            E.badinput('If you wish to filter by time, the input data must represent 366 days');
                        end
                    end
                    
                    %gridded_data = nanmean(gridded_data,3);
                    gridded_data = nanmedian(gridded_data,3);
                    gridded_data = gridded_data(xx_centers,yy_centers);
                end
                
                
                
                sat_lon = geos_chem_centers('2x25');
                sat_lon = sat_lon(xx_centers);
        end
        
        % Make the meridional average
        %sat_avg = nanmean(gridded_data,2)';
        sat_avg = nanmedian(gridded_data,2)';
        
    end

    function gridded_daily_omno2d = grid_omno2d_daily(gc_loncorn, gc_latcorn)
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/OMI Utils');
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        omno2d_path = '/Volumes/share-sat/SAT/OMI/OMNO2d';
        if ~exist(omno2d_path,'dir')
            E.dir_dne('omno2d_path')
        end
        [omno2d_lon, omno2d_lat] = omno2d_centers;
        
        fname = sprintf('OMI-Aura_L3-OMNO2d_%04dm%02d%02d_v003*.he5', year(dates(d)), month(dates(d)), day(dates(d)));
        F = dir(fullfile(omno2d_path, sprintf('%04d',year(dates(d))),fname));
        if isempty(F) && ~ignore_missing
            E.filenotfound(fname);
        elseif numel(F) > 1
            E.toomanyfiles(fname);
        end
        
        gridded_daily_omno2d = nan([size(gc_loncorn)-1, 366]);
        
        for d=1:366
            sdate = modis_day_to_date(d,2012);
            
            hinfo = h5info(fullfile(omno2d_path, sprintf('%04d',year(sdate)), F(1).name));
            todays_no2 = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2TropCloudScreened'));
            todays_weight = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'Weight'));
            
            % Ignore negative values - they are either fill values or unphysical
            todays_weight(todays_no2<0) = 0;
            
            gridded_daily_omno2d(:,:,d) = grid_omno2d_to_gc(todays_no2, omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn, todays_weight);
        end
    end

    function [gridded_monthly_omno2d, monthly_omno2d] = grid_omno2d_monthly(gc_loncorn, gc_latcorn)
        % Does monthly time averages of OMNO2d data then grids it to GC
        % resolution.
        [omno2d_lon, omno2d_lat] = omno2d_centers;
        monthly_omno2d = nan(size(omno2d_lon));
        gridded_monthly_omno2d = nan(size(gc_loncorn)-1);
        for m=1:12
            sdate = sprintf('2012-%02d-01',m);
            edate = sprintf('2012-%02d-%02d',m,eomday(2012,m));
            monthly_omno2d(:,:,m) = omno2d_timeavg(sdate,edate);
            gridded_monthly_omno2d(:,:,m) = grid_omno2d_to_gc(monthly_omno2d(:,:,m), omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn);
        end
        
    end

    function gridded_omno2d = grid_omno2d_to_gc(omno2d_no2, omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn, omno2d_weights)
        % Averages OMNO2d data to GEOS-Chem resolution. Needs OMNO2d data,
        % longitude, and latitude (all should be the same size) and the
        % GEOS-Chem corner points to average to.
        
        % Check input
        if ndims(omno2d_no2) ~= ndims(omno2d_lon) || ndims(omno2d_no2) ~= ndims(omno2d_lat) || any(size(omno2d_no2) ~= size(omno2d_lon)) || any(size(omno2d_no2) ~= size(omno2d_lat))
            E.sizeMismatch('omno2d_lon','omno2d_lat','omno2d_no2')
        end
        if ~exist(omno2d_weights)
            omno2d_weights = ones(size(omno2d_no2));
        end
        
        gridded_omno2d = nan(size(gc_loncorn)-1); % one smaller to switch from # of corners to # of grid cells
        
        % Loop over GC grid cells and find which OMNO2d data belongs in it
        for a=1:size(gridded_omno2d,1)
            if DEBUG_LEVEL > 1; fprintf('%d ',a); end
            for b=1:size(gridded_omno2d,2)
                xx = omno2d_lon >= gc_loncorn(a,b) & omno2d_lon < gc_loncorn(a+1,b+1) & omno2d_lat >= gc_latcorn(a,b) & omno2d_lat < gc_latcorn(a+1,b+1);
                gridded_omno2d(a,b) = nansum(omno2d_no2(xx) .* omno2d_weights(xx)) / nansum(omno2d_weights(xx));
            end
        end
        if DEBUG_LEVEL > 1; fprintf('\n'); end
    end

    function [gridded_mls, count] = grid_daily_mls_data_to_gc(gc_loncorn, gc_latcorn, specie, omi_overpass)
        % Calculates a running average of MLS data for the given GC
        % longitude and latitude corners.  Intended to be called from
        % sat_meridional_avg() or externally to pre-grid data and save
        % yourself time. (In that case, pass a full gc_loncorn/latcorn
        % matrix - don't cut it down into regions.)
        %
        % One optional argument, set to true to limit to approximately OMI
        % overpass time (will restrict to 11:00 - 15:00 solar time).
        
        if ~exist('omi_overpass','var')
            omi_overpass = false;
        end
        
        switch lower(specie)
            case 'o3'
                filepath = '/Volumes/share-sat/SAT/MLS/O3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [350 200];
            case 'hno3'
                filepath = '/Volumes/share-sat/SAT/MLS/HNO3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                %pressure_range = [250 200];
                %pressure_range = [250 140];
                pressure_range = [200 140];
        end
        

        gridded_mls = nan([size(gc_loncorn)-1, 366]);
        count = zeros([size(gc_loncorn)-1, 366]);

        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        
        fnames = {files.name};
        
        % Do the gridding and averaging
        d_ind = 0;
        for d=datenum('2012-01-01'):datenum('2012-12-31')
            if DEBUG_LEVEL > 1
                fprintf('Binning %s\n',datestr(d));
            end
            d_ind = d_ind + 1;
            fileglob = sprintf('MLS-Aura_L2GP-%s_v04-2[0-9]-c01_2012d%03d.he5',upper(specie),modis_date_to_day(d));
            filename = glob(fnames, fileglob);
            if isempty(filename)
                continue
            end
            Data = read_mls_data(fullfile(filepath,filename{1}));
            lon = Data.Longitude;
            lat = Data.Latitude;
            solar_time = Data.LocalSolarTime;
            
            % Find which levels of the MLS retrieval are in the desired
            % pressure range
            pp = Data.Pressure >= min(pressure_range) & Data.Pressure <= max(pressure_range);
            
            % Restrict (more or less) to OMI overpass time based on the
            % local solar time of the measurement. Mainly intended to
            % remove nighttime rather than completely limit to OMI
            % overpass.
            if omi_overpass
                tt = solar_time >= 11 & solar_time <= 15;
                Data.(upper(specie))(:,~tt) = nan;
            end
            
            for a = 1:size(gridded_mls,1)
                for b = 1:size(gridded_mls,2)
                    % Find the elements from the MLS data that fall in each
                    % GC grid cell and add them to the running average
                    xx = lon >= gc_loncorn(a,b) & lon < gc_loncorn(a+1,b+1) & lat >= gc_latcorn(a,b) & lat < gc_latcorn(a+1,b+1);
                    if sum(xx) < 1
                        continue
                    end
                    c = sum(xx);
                    data = nansum2(Data.(upper(specie))(pp,xx));
                    
                    gridded_mls(a, b, d_ind) = nansum2([gridded_mls(a,b) * count(a,b), data]) / (count(a,b) + c);
                    count(a, b, d_ind) = count(a,b) + c;
                end
            end
        end
    end

    function gridded_mls = grid_mls_data_to_gc(gc_loncorn, gc_latcorn, specie)
        % Calculates a running average of MLS data for the given GC
        % longitude and latitude corners.  Intended to be called from
        % sat_meridional_avg() or externally to pre-grid data and save
        % yourself time. (In that case, pass a full gc_loncorn/latcorn
        % matrix - don't cut it down into regions.)
        
        switch lower(specie)
            case 'o3'
                filepath = '/Volumes/share-sat/SAT/MLS/O3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [350 200];
            case 'hno3'
                filepath = '/Volumes/share-sat/SAT/MLS/HNO3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [250 200];
        end
        
        if temporal_avg
            gridded_mls = nan(size(gc_loncorn)-1);
            count = zeros(size(gc_loncorn)-1);
        else
            gridded_mls = nan([size(gc_loncorn)-1, numel(files)]);
            count = zeros([size(gc_loncorn)-1, numel(files)]);
        end
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        
        % Do the gridding and averaging
        for f=1:numel(files)
            if DEBUG_LEVEL > 1
                fprintf('Binning file %d of %d\n',f,numel(files));
            end
            Data = read_mls_data(fullfile(filepath,files(f).name));
            lon = Data.Longitude;
            lat = Data.Latitude;
            
            % Find which levels of the MLS retrieval are in the desired
            % pressure range
            pp = Data.Pressure >= min(pressure_range) & Data.Pressure <= max(pressure_range);
            
            for a = 1:size(gridded_mls,1)
                for b = 1:size(gridded_mls,2)
                    % Find the elements from the MLS data that fall in each
                    % GC grid cell and add them to the running average
                    xx = lon >= gc_loncorn(a,b) & lon < gc_loncorn(a+1,b+1) & lat >= gc_latcorn(a,b) & lat < gc_latcorn(a+1,b+1);
                    if sum(xx) < 1
                        continue
                    end
                    c = sum(xx);
                    d = nansum(Data.(upper(specie))(pp,xx));
                    
                    if temporal_avg
                        gridded_mls(a,b) = nansum2([gridded_mls(a,b) * count(a,b), d]) / (count(a,b) + c);
                        count(a,b) = count(a,b) + c;
                    else
                        gridded_mls(a,b,f) = nansum2([gridded_mls(a,b) * count(a,b), d]) / (count(a,b) + c);
                        count(a,b,f) = count(a,b) + c;
                    end
                end
            end
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED PLOTTING FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_dcolumn_percentiles(case1sat, case1crit, case2sat, case2crit, struct_ind, varargin)
        % Plots a timeseries of the difference in total columns between two
        % cases. Requires structures with the information about the columns
        % as the first and third arguments, and logical matrices indicating
        % where lightning NOx is >60% of the emission as the second and
        % fourth. The next to last argument is the index for which field in
        % the structures to plot. The last argument is a boolean
        % determining whether to plot the absolute or relative difference
        % in columns. This is set by the internal function call and does
        % not need to be passed as a user argument to the main function.
        % The difference is defined as first - second or (first -
        % second)/second.
        
        % Parse varargs. If one is the string 'relative', the next one is
        % the relative boolean. There should be 0 or 2 left after that, if
        % 2, then those are the lat and lon limits.
        [relative, ~, lonlim, latlim] = parse_varargs(varargin);
        
        if ~ismember('Columns',fieldnames(case1sat)) %isfield wasn't behaving
            case1sat = integrate_geoschem_profile(case1sat,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2sat))
            case2sat = integrate_geoschem_profile(case2sat,1e-9);
        end
        
        [glon, glat] = geos_chem_centers(size(case1sat(1).dataBlock));
        xx = glon >= min(lonlim) & glon <= max(lonlim);
        yy = glat >= min(latlim) & glat <= max(latlim);
        
        case1col = case1sat(struct_ind).Columns;
        case2col = case2sat(struct_ind).Columns;
        
        case1col(~case1crit) = nan;
        case2col(~case2crit) = nan;
        
        case1col = case1col(xx,yy,:);
        case2col = case2col(xx,yy,:);
        if ~relative
            dcol = case1col - case2col;
        else
            dcol = (case1col - case2col) ./ case2col * 100;
        end
        
        plot_percentile_timeseries(dcol, case1sat(struct_ind).tVec);
        if ~relative
            ylabel('Absolute column differences (molec./cm^2)');
        else
            ylabel('Percent difference in column')
        end
    end

    function plot_level_dconc(case1sat, case1crit, case2sat, case2crit, struct_ind, plevels, varargin)
        % Most inputs the same as plot_dcolumn_percentiles, plevels is a
        % vector of level indicies to make plots for, 1 per level. 
        
        [relative, ndens_bool, lonlim, latlim] = parse_varargs(varargin);
               
        [glon, glat] = geos_chem_centers(size(case1sat(1).dataBlock));
        xx = glon >= min(lonlim) & glon <= max(lonlim);
        yy = glat >= min(latlim) & glat <= max(latlim);
        
        
        titlestr = input('Give a title for these plots with a %s for the pressure level: ', 's');
        caseinds = get_data_inds(case1sat);
        for a=1:numel(plevels)
            case1conc = squeeze(case1sat(struct_ind).dataBlock(xx,yy,plevels(a),2:end))*1e-9;
            case1ndens = squeeze(case1sat(caseinds.ndens).dataBlock(xx,yy,plevels(a),2:end))*1e-6;
            case2conc = squeeze(case2sat(struct_ind).dataBlock(xx,yy,plevels(a),2:end))*1e-9;
            case2ndens = squeeze(case2sat(caseinds.ndens).dataBlock(xx,yy,plevels(a),2:end))*1e-6;
            if ndens_bool
                case1conc = case1conc;% .* case1ndens;
                case2conc = case2conc;% .* case2ndens;
            end
            case1conc(~case1crit(xx,yy,2:end)) = nan;
            case2conc(~case2crit(xx,yy,2:end)) = nan;

            
            if ~relative
                dconc = case1conc - case2conc;
            else
                dconc = (case1conc - case2conc)./case2conc * 100;
            end
            
            % We don't want to include any results from above the
            % tropopause, so set those to NaN.
            tpbool = floor(case1sat(caseinds.tplev).dataBlock(xx,yy,2:end)) <= plevels(a);
            dconc(tpbool) = nan;
            
            plot_percentile_timeseries(dconc, case1sat(struct_ind).tVec(2:end));
            
            if ~relative && ndens_bool
                ylabel('Absolute concentration differences (molec./cm^3)');
            elseif ~relative && ~ndens_bool
                ylabel('Absolute VMR differences (ppp)');
            elseif relative && ndens_bool
                ylabel('Percent difference in num. density')
            elseif relative && ~ndens_bool
                ylabel('Percent difference in VMR');
            end
            
            titlestr = tex_in_printf(titlestr, 2);
            titlestr2 = sprintf(titlestr, '(trop only, Plev = %d, P range=[%.1f, %.1f])');
            sz = size(case1sat(caseinds.psurf).dataBlock(xx,yy,:,:));
            % List the range of pressures found for a given vertical
            % coordinate. However, we need to skip the first day because
            % any ND51 diagnostics that output not at midnight give weird
            % results on the first day
            minpres = min(reshape(case1sat(caseinds.psurf).dataBlock(xx,yy,plevels(a),2:end),sz(1)*sz(2)*(sz(4)-1),1));
            maxpres = max(reshape(case1sat(caseinds.psurf).dataBlock(xx,yy,plevels(a),2:end),sz(1)*sz(2)*(sz(4)-1),1));
            title(sprintf(titlestr2, plevels(a), minpres, maxpres));
        end
    end

    function plot_lnox_col_enhnc_diff(case1lnox, case1nolnox, case1crit, case2lnox, case2nolnox, case2crit, struct_ind, varargin)
        % Plots quantile differences in lightning enhancement over time.
        % case[12][lnox|nolnox] should be the same structures as
        % case[12]sat in the previous two functions, just one with and one
        % without lightning. case[12]crit likewise is the same logical
        % matrix of profiles with significant lightning emissions (e.g. 60%
        % of NO emissions from lightning, BB, anthro due to lightning).
        % struct_ind and relative are the same as before.
        
        [relative, ndens_bool, lonlim, latlim] = parse_varargs(varargin);
        
        if ~ismember('Columns',fieldnames(case1lnox)) %isfield wasn't behaving
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case1lnox\n'); end
            case1lnox = integrate_geoschem_profile(case1lnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case1nolnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case1nolnox\n'); end
            case1nolnox = integrate_geoschem_profile(case1nolnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2lnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case2lnox\n'); end
            case2lnox = integrate_geoschem_profile(case2lnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2nolnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case2nolnox\n'); end
            case2nolnox = integrate_geoschem_profile(case2nolnox,1e-9);
        end
        
        case1lnox_col = case1lnox(struct_ind).Columns;
        case1nolnox_col = case1nolnox(struct_ind).Columns;
        case2lnox_col = case2lnox(struct_ind).Columns;
        case2nolnox_col = case2nolnox(struct_ind).Columns;
        
        case1lnox_col(~case1crit) = nan;
        case1nolnox_col(~case1crit) = nan;
        case2lnox_col(~case2crit) = nan;
        case2nolnox_col(~case2crit) = nan;
        
        dcol = (case1lnox_col - case1nolnox_col) - (case2lnox_col - case2nolnox_col);
        if relative
            dcol = dcol ./ (case2lnox_col - case2nolnox_col) * 100;
        end
        
        plot_percentile_timeseries(dcol, case1lnox(struct_ind).tVec);
        if ~relative
            ylabel('Absolute enhancement differences (molec./cm^2)');
        else
            ylabel('Percent enhancement difference in column')
        end
    end

    function plot_percentile_timeseries(mat_in, tVec)
        % Assume the third dimension is time in the column matrices
        sz = size(mat_in);
        casequant = zeros(sz(3),5);
        
        for a=1:sz(3)
            sub = mat_in(:,:,a);
            casequant(a,:) = quantile(sub(:),[0.05, 0.25, 0.5, 0.75, 0.95]);
        end
        
        figure;
        plot(tVec, casequant);
        datetick('x');
        legend('0.05','0.25','0.50','0.75','0.95');
        set(gca,'ygrid','on');
    end

    function [avg_data, glon, region, p_max, column_bool, fig] = meridional_avg(dataStruct, criteria, varargin)

        
        % Parse the optional arguments
        fig = [];
        column_bool = false;
        timelim_bool = false;
        p_max = [];
        p_mat = [];
        region = '';
        tp_mask = [];
        for a=1:numel(varargin)
            if isnumeric(varargin{a}) && (isscalar(varargin{a}) || (isvector(varargin{a}) && numel(varargin{a}) == 2))
                p_max = varargin{a};
            elseif isstruct(varargin{a}) && ~isempty(regexpi(varargin{a}.fullName, 'PSURF', 'once'))
                p_mat = varargin{a}.dataBlock;
            elseif (isstruct(varargin{a}) && strcmpi(varargin{a}.fullName, 'Tropopause level')) || islogical(varargin{a})
                tp_mask = varargin{a};
            elseif ischar(varargin{a}) && strcmpi('columns',varargin{a})
                column_bool = true;
            elseif ~ischar(varargin{a}) && ishandle(varargin{a}) && strcmpi(get(varargin{a},'type'),'figure')
                fig = varargin{a};
            elseif ischar(varargin{a}) && ismember(varargin{a}, {'na','sa','naf','saf','seas','neur'})
                if ~isempty(region)
                    E.badinput('Multiple regions specified ("%s" and "%s"). Only specify one.', region, varargin{a});
                end
                region = varargin{a};
            elseif ischar(varargin{a}) && strcmpi(varargin{a}, 'timelim')
                timelim_bool = true;
            end
        end
        
        [lonlim, latlim, timemask] = define_regions(region);
        
        % Make sure the input structure is correct (has everything needed)
        req_fields = {'dataUnit','fullName'};
        if column_bool
            req_fields{end+1} = 'Columns';
        else
            req_fields{end+1} = 'dataBlock';
        end
        
        if ~isstruct(dataStruct) || ~isscalar(dataStruct)
            E.badinput('dataStruct must be a scalar structure')
        end
        flds = ismember(req_fields, fieldnames(dataStruct));
        if any(~flds)
            E.badinput('dataStruct must have the fields %s', strjoin(req_fields(~flds),', '));
        end
        
        % Read in Columns or dataBlock. If a level isn't specified and
        % we're reading in dataBlock, the user forgot.
        if column_bool
            data = dataStruct.Columns;
        else
            if isempty(p_max)
                E.badinput('Max pressure must be specified if plotting a concentration')
            elseif isempty(p_mat)
                E.badinput('A matrix of pressures must be give if plotting concentrations')
            elseif isempty(tp_mask) && isscalar(p_max)
                E.badinput('A troposphere mask must be given if plotting a concentration and only a max pressure given')
            elseif isstruct(tp_mask) && ~strcmpi(tp_mask.fullName,'Tropopause level')
                E.badinput('The tropopause mask is a structure but does not seem to refer to the tropopause level')
            end
            
            % Set stratospheric boxes to be NaN
            if isscalar(p_max)
                if isstruct(tp_mask)
                    tp_mask = trop_mask(size(dataStruct.dataBlock), tp_mask);
                end
                dataStruct.dataBlock(~tp_mask) = nan;
            end

            % Average the concentrations over the requested levels. Set
            % boxes below the max pressure to be NaN.
            if ndims(dataStruct.dataBlock) ~= ndims(p_mat) || any(size(dataStruct.dataBlock) ~= size(p_mat))
                E.badinput('Pressure and concentration matrices must be the same size')
            end
            
            if isscalar(p_max)
                pp = p_mat <= p_max;
            else
                pp = p_mat <= max(p_max) & p_mat >= min(p_max);
            end
            dataStruct.dataBlock(~pp) = nan;
            data = squeeze(nanmean(dataStruct.dataBlock,3));
        end
        
        % criteria must be a logical matrix with the same lat/lon/time
        % dimensions as dataBlock or Columns in the input structure.
        if any(size(data) ~= size(criteria))
            E.badinput('Size of input data and criteria do not match')
        end
        
        % Set times outside the time range to NaN so they are not included
        % in the average. Do this only if the 'timelim' argument has been
        % passed and if "data" is 366 elements long in the third dimension,
        % as the time limits assume so (i.e. 366 days of 2012).
        if timelim_bool
            if size(data,3) == 366
                data(:,:,~timemask) = nan;
            else
                E.badinput('To filter by time, the input matrix must represent 366 days')
            end
        end
        
        % Needed for x ticks and to limit to the region
        [glon, glat] = geos_chem_centers(size(data));
        
        % Average the data meridionally (i.e. along latitude). Ignore boxes
        % where the criterion is false. Martin seems to time average first,
        % then average latitudinally, so we'll do the same. Remove boxes
        % outside the area of interest first.
        
        data(~criteria) = nan;
        
        xx = glon >= lonlim(1) & glon <= lonlim(2);
        yy = glat >= latlim(1) & glat <= latlim(2);
        data = data(xx,yy,:);
        
        %avg_data = nanmean(nanmean(data,3),2);
        avg_data = nanmedian(nanmedian(data,3),2);
        glon = glon(xx);
    end

    function plot_meridional_avg(dataStruct, criteria, varargin)
        % Will produce a figure like Fig. 5 in Randall Martin's 2007
        % lightning NOx paper, where the NO2 columns or O3/HNO3
        % concentration are averaged meridionally for places that fit the
        % criteria. Req. 2 arguments: dataStruct must be one structure
        % output from read_geos_output or cleaned up from the ts_satellite
        % python reading code, and the criteria matrix (a logical matrix
        % true wherever the data should be included, i.e lightning >60% of
        % emissions). A troposphere mask - this can either be a logical
        % matrix the same size as the dataBlock or the tropopause level
        % output structure from read_geos_output - is required if doing
        % levels instead of columns.
        %
        % If plotting a level, you'll need to include a matrix of pressures
        % and the maximum pressure to plot.
        %
        % Optional arguments: a handle to a figure will plot on that figure
        % instead of creating a new one. The string 'columns' will cause it
        % to average columns rather than concentrations. Other strings will
        % specify what lon/lat limits to apply: 'na' (north america), 'sa'
        % (south america), 'naf' (north africa), 'saf' (south africa), and
        % 'seas' (southeast asia). The string 'timelim' will limit the
        % region to its three summer months (JJA for N. Hemisphere or DJF
        % for SH).
        
        [avg_data, glon, region, p_max, column_bool, fig] = meridional_avg(dataStruct, criteria, varargin{:});
        
        % Create or select the figure as necessary
        if isempty(fig)
            figure;
        else
            figure(fig);
            hold on
        end
        
        plot(glon', avg_data);
        
        if isempty(fig)
            species_name = regexprep(dataStruct.fullName,'[_^]',' ');
            xlabel('Longitude')
            if column_bool
                ylabel(sprintf('%s columns (molec. cm^{-2})', species_name));
            else
                ylabel(sprintf('%s (%s)', species_name, regexprep(dataStruct.dataUnit, '[_^]', ' ')));
            end
            if ~isempty(region)
                region_title = sprintf(' over %s',upper(region));
            else
                region_title = '';
            end
            if ~isempty(p_max)
                level_title = sprintf(' for P \\\\leq %d', p_max);
            else
                level_title = '';
            end
            species_name = strsplit(species_name);
            title(sprintf('%s Meridonal avg.%s%s', species_name{2}, region_title, level_title));
            set(gca,'fontsize',14);
        end
    end

    function add_sat_mer_avg_to_fig(fig, specie, region, gridded_data, lnox_crit, varargin)
        % Will call sat_meridional_avg to calculate the meridional average
        % of satellite data for the given region then add that to the
        % figure given by the handle fig.  The figure should have a legend
        % already.  If you have pregridded data, pass as the final argument
        % to save time gridding.  This WILL NOT check that you are adding
        % a like species to the figure, so it will happily say put NO2
        % columns on an HNO3 concentration plot or put data from N. Am. on
        % SE Asia. Just warning you...
        %
        % Some optional arguments: the string 'noregrid' will keep the data
        % in its native gridding, and the string 'timelim' will limit the
        % data for each region to that region's three summer months (JJA or
        % DJF).
        
        % Check & parse input
        if ~ishandle(fig) || ~strcmpi(fig.Type, 'figure')
            E.badinput('fig must be a handle to a figure')
        end
        if ~exist('gridded_data','var')
            gridded_data = []; % placeholder, sat_meridional_avg will recognize this as indicating gridding must be done.
        end
        if ~exist('lnox_crit','var')
            lnox_crit = [];
        end
        noregrid_bool = false;
        timelim_bool = false;
        if numel(varargin) > 0
            for a = 1:numel(varargin)
                if ischar(varargin{a}) && strcmpi(varargin{a}, 'noregrid')
                    noregrid_bool = true;
                elseif ischar(varargin{a}) && strcmpi(varargin{a}, 'timelim')
                    timelim_bool = true;
                end
            end
        end
        
        ax = findobj(fig,'Type','axes');
        leg = findobj(fig,'Type','legend');
        
        [sat_lon, sat_avg] = sat_meridional_avg(specie, region, gridded_data, lnox_crit, noregrid_bool, timelim_bool);
        
        line(sat_lon, sat_avg, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--', 'parent', ax);
        %lstrings = leg.String;
        %lstrings{end+1} = 'Satellite';
        %legend(lstrings{:});
    end

    function [ avg_table ] = avg_meridional_avg_table(dataStruct_new, criteria_new, dataStruct_old, criteria_old, varargin)
       % Takes the meridional averages for all 6 regions and then averages
       % over the longitudes, giving a table of absolute and relative
       % average differences in meridional averages by region. See
       % plot_meridional_avg for a description of the inputs, with two
       % differences: 
       %    One, you must input two data structures (the "new" and "old"
       %    cases) as this is computing a difference.
       %
       %    Two, do not input a region as this function will reject it
       %    before calling meridional_avg, as this function will supply the
       %    regions itself.
       
       % This will be used in the table, so it should be a column vector
       regions = {'na';'sa';'naf';'saf';'seas';'neur'};
       
       for a=1:numel(varargin)
           if ischar(varargin{a})
               xx = ismember(varargin{a}, regions);
               
               if any(xx)
                   E.badinput('A region ("%s") was detected in input; remove this and try again. This function will iterate over all regions internally', varargin{a});
               end
           end
       end
       
       AbsoluteDifference = nan(size(regions));
       RelativeDifference = nan(size(regions));
       
       for a=1:numel(regions)
           avg_data_new = meridional_avg(dataStruct_new, criteria_new, varargin{:}, regions{a});
           avg_data_old = meridional_avg(dataStruct_old, criteria_old, varargin{:}, regions{a});
           
           AbsoluteDifference(a) = nanmean(avg_data_new - avg_data_old);
           RelativeDifference(a) = nanmean((avg_data_new - avg_data_old) ./ avg_data_old * 100);
       end
       
       avg_table = table(AbsoluteDifference, RelativeDifference, 'RowNames', regions);
    end

    function [omno2_no2, omno2_std, domino_no2, domino_std, omno2_gc, domino_gc] = compare_regions_to_sat(varargin)
        % Will generate the (hopefully) final figure for Ben and my paper.
        % The current idea I have is to make two panels, one for OMNO2 one
        % for DOMINO that compare the average satellite columns for the
        % four "good" satellite regions (S.Am., S.Af., N.Af., and SE Asia)
        % to the AK-corrected base, HendwMPN-PNA-N2O5 and HendwMPN-PNA-N2O5
        % +33% cases. There will need to be two plots I think because the
        % AKs are different for each product. So each region will have 8
        % points, four per panel.
        
        if ismember('sdom',varargin)
            sdom = true;
            fprintf('Using std. dev. of mean\n');
        else
            sdom = false;
            fprintf('Using just std. dev.\n');
        end

        % Load the production based masks that we'll need to cut down the
        % GC and satellite matices to the lightning heavy observations. We
        % will use the lnox > 60% from the base runs for all runs so that
        % we are comparing the same grid cells in every case.
        M = load(fullfile('/Users/Josh/Documents/MATLAB/','MPN Project','Workspaces','Pickering Parameterization','Daily24hrs','All2012Daily','HendwMPN_PNA_N2O5-Pickp33-ProdFrac.mat'), 'hendwmpn_pna_n2o5_pickp33_lnoxgt60','anthro_lt15e10');
        lnox_mask = M.hendwmpn_pna_n2o5_pickp33_lnoxgt60;
        anthro_mask = M.anthro_lt15e10;
        
        % Then load the satellite data and average it down for each
        % regions.
        
        omno2_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/OMNO2/2.5x2.0-avg';
        domino_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/DOMINO/2.5x2.0-avg';
        
        fprintf('Loading OMNO2 data... Please be patient...   ');
        F_omno2 = dir(fullfile(omno2_path,'OMI*.mat'));
        omno2_data = nan(144,91,numel(F_omno2));
        omno2_weight = nan(144,91,numel(F_omno2));
        for a=1:numel(F_omno2)
            O = load(fullfile(omno2_path, F_omno2(a).name));
            if a == 1
                omno2_lon = O.GC_avg.Longitude;
                omno2_lat = O.GC_avg.Latitude;
            end
            omno2_data(:,:,a) = O.GC_avg.ColumnAmountNO2Trop;
            omno2_weight(:,:,a) = O.GC_avg.Weight;
            clear('O')
        end
        fprintf('Done.\n');
        omno2_mask = ~isnan(omno2_data);
        omno2_data(~lnox_mask | ~anthro_mask) = nan;
        
        fprintf('Loading DOMINO data... Please be patient...   ');
        F_domino = dir(fullfile(domino_path, 'OMI*.mat'));
        domino_data = nan(144,91,numel(F_domino));
        domino_weight = nan(144,91,numel(F_domino));
        for a=1:numel(F_domino)
            D = load(fullfile(domino_path, F_domino(a).name));
            if a == 1
                domino_lon = D.GC_avg.Longitude;
                domino_lat = D.GC_avg.Latitude;
            end
            domino_data(:,:,a) = D.GC_avg.TroposphericVerticalColumn;
            domino_weight(:,:,a) = D.GC_avg.Weight;
            clear('D')
        end
        domino_mask = ~isnan(domino_data);
        domino_data(~lnox_mask | ~anthro_mask) = nan;
        fprintf('Done.\n');
        
        fprintf('Dividing sat data into regions...   ');
        regions = {'sa','saf','naf','seas'};
        omno2_no2 = make_empty_struct_from_cell(regions);
        omno2_std = make_empty_struct_from_cell(regions);
        domino_no2 = make_empty_struct_from_cell(regions);
        domino_std = make_empty_struct_from_cell(regions);
        for r=1:numel(regions)
            [xx,yy] = region_xxyy(regions{r});
            [~,~,timeind] = define_regions(regions{r});
            omno2_no2.(regions{r}) = nanmean(reshape(omno2_data(xx,yy,timeind),1,[]));
            omno2_std.(regions{r}) = nanstd(reshape(omno2_data(xx, yy, timeind),1,[]));
            
            domino_no2.(regions{r}) = nanmean(reshape(domino_data(xx, yy, timeind),1,[]));
            domino_std.(regions{r}) = nanstd(reshape(domino_data(xx, yy, timeind),1,[]));
        end
        fprintf('Done.\n');
        
        % Now load each of the 3 cases from GEOS-Chem for each of the
        % products' AKs and get the same regions.
        gc_path = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/DailyOMI/AllAKsDaily';
        chems = {'JPLnoMPN','HendwMPN-PNA-N2O5','HendwMPN-PNA-N2O5'};
        omno2_ak_files = {'JPLnoMPN-Pickp0-OMI_NO2_aks.mat','HendwMPN_PNA_N2O5-Pickp0-OMI_NO2_OMI_aks.mat','HendwMPN_PNA_N2O5-Pickp33-OMI_NO2_OMI_aks.mat'};
        dom_ak_files = {'JPLnoMPN-Pickp0-OMI_NO2_DOMINO_aks.mat','HendwMPN-PNA-N2O5-Pickp0-OMI_NO2_DOMINO_aks.mat','HendwMPN-PNA-N2O5-Pickp33-OMI_NO2_DOMINO_aks.mat'};
        fields = {'Base','Final','Final33'};
        gc_data = make_empty_struct_from_cell(fields);
        omno2_gc = make_empty_struct_from_cell(regions,gc_data);
        domino_gc = make_empty_struct_from_cell(regions,gc_data);
        for a=1:numel(omno2_ak_files)
            fprintf('Loading %s...   ', omno2_ak_files{a});
            O = load(fullfile(gc_path, chems{a}, omno2_ak_files{a}));
            fns = fieldnames(O);
            if numel(fns) > 1; warning('Multiple variables in %s, reading %s only',omno2_ak_files{a},fns{1}); end
            gc = O.(fns{1}).Columns;
            gc(~lnox_mask | ~anthro_mask | ~omno2_mask) = nan;
            for r=1:numel(regions)
                [xx,yy] = region_xxyy(regions{r});
                [~,~,timeind] = define_regions(regions{r});
                omno2_gc.(regions{r}).(fields{a}).mean = nanmean(reshape(gc(xx,yy,timeind),1,[]));
                omno2_gc.(regions{r}).(fields{a}).std = nanstd(reshape(gc(xx,yy,timeind),1,[]));
            end
            fprintf('Done.\n');
            clear('O','gc')
            
            fprintf('Loading %s...   ', dom_ak_files{a});
            D = load(fullfile(gc_path, chems{a}, dom_ak_files{a}));
            fns = fieldnames(D);
            if numel(fns) > 1; warning('Multiple variables in %s, reading %s only',dom_ak_files{a},fns{1}); end
            gc = D.(fns{1}).Columns;
            gc(~lnox_mask | ~anthro_mask | ~domino_mask) = nan;
            for r=1:numel(regions)
                [xx,yy] = region_xxyy(regions{r});
                [~,~,timeind] = define_regions(regions{r});
                domino_gc.(regions{r}).(fields{a}).mean = nanmean(reshape(gc(xx,yy,timeind),1,[]));
                domino_gc.(regions{r}).(fields{a}).std = nanstd(reshape(gc(xx,yy,timeind),1,[]));
            end
            fprintf('Done.\n')
            clear('D','gc')
        end
        
        if sdom
            sqrtn_o = sqrt(sum(lnox_mask(:) & anthro_mask(:) & omno2_mask(:)));
            sqrtn_d = sqrt(sum(lnox_mask(:) & anthro_mask(:) & domino_mask(:)));
        else
            sqrtn_o = 1;
            sqrtn_d = 1;
        end

        % Finally make the figures. Two panels: OMNO2 on top, DOMINO on
        % bottom.
        
        figure; 
        subplot(2,1,1)
        
        cols = {'r','b',[0 0.5 0]};
        marks = {'x','^','s'};
        
        for r=1:numel(regions)
            x = (r-1)*2+1;
            l(1) = line(x-0.5,omno2_no2.(regions{r}),'color','k','marker','o','linewidth',2);
            scatter_errorbars(x-0.5, omno2_no2.(regions{r}), omno2_std.(regions{r})/sqrtn_o, 'color','k','linewidth',2);
            for a=1:numel(fields)
                dx = a*0.25 - 0.5;
                l(a+1) = line(x+dx, omno2_gc.(regions{r}).(fields{a}).mean, 'color', cols{a}, 'marker', marks{a}, 'linewidth', 2);
                scatter_errorbars(x+dx, omno2_gc.(regions{r}).(fields{a}).mean, omno2_gc.(regions{r}).(fields{a}).std/sqrtn_o, 'color', cols{a}, 'linewidth',2);
            end
        end
        set(gca,'xtick',[1 3 5 7]);
        set(gca,'xticklabels',{'S. Am.','S. Af.','N. Af.', 'SE Asia'});
        set(gca,'fontsize',16);
        legend(l',{'OMNO2','Base case','Final case','Final case +33%'});
        
        subplot(2,1,2);
        for r=1:numel(regions)
            x = (r-1)*2+1;
            l(1) = line(x-0.5,domino_no2.(regions{r}),'color','k','marker','o','linewidth',2);
            scatter_errorbars(x-0.5, domino_no2.(regions{r}), domino_std.(regions{r})/sqrtn_d, 'color','k','linewidth',2);
            for a=1:numel(fields)
                dx = a*0.25 - 0.5;
                l(a+1) = line(x+dx, domino_gc.(regions{r}).(fields{a}).mean, 'color', cols{a}, 'marker', marks{a}, 'linewidth', 2);
                scatter_errorbars(x+dx, domino_gc.(regions{r}).(fields{a}).mean, domino_gc.(regions{r}).(fields{a}).std/sqrtn_d, 'color', cols{a}, 'linewidth',2);
            end
        end
        set(gca,'xtick',[1 3 5 7]);
        set(gca,'xticklabels',{'S. Am.','S. Af.','N. Af.', 'SE Asia'});
        set(gca,'fontsize',16);
        legend(l',{'DOMINO','Base case','Final case','Final case +33%'});
    end
end

