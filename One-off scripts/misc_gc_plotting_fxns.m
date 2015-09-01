function [ varargout ] = misc_gc_plotting_fxns( plttype, varargin )
%GC_PLOTTING_FXNS Misc. plots to make for GC outputs
%   Rather than trying to do this all in the command window, or having a
%   zillion different .m files all over the place, I'm just going to
%   collect them all here.
%
%   Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2015

DEBUG_LEVEL = 1;
E = JLLErrors;

plttype = lower(plttype);
switch plttype
    case 'subfrac_lnox'
        [varargout{1}, varargout{2}] = subfrac_lnox(varargin{1});
    case 'trop_mask'
        [varargout{1}] = trop_mask(varargin{1:2});
    case 'region_xxyy'
        [varargout{1}, varargout{2}] = region_xxyy(varargin{1});
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

    function [relative, ndens_bool, lonlim, latlim] = parse_varargs(vargs)
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

    function [lonlim, latlim] = define_regions(region)
        switch region
            case 'na'
                lonlim = [-120, -65];
                latlim = [20, 60];
            case 'sa'
                lonlim = [-77, -39];
                latlim = [-35, 10];
            case 'naf'
                lonlim = [-15, 48];
                latlim = [3, 25];
            case 'saf'
                lonlim = [10, 48];
                latlim = [-30, 3];
            case 'seas'
                lonlim = [95, 146];
                latlim = [-9, 26];
            case 'neur'
                lonlim = [60, 130];
                latlim = [30, 68];
            otherwise
                lonlim = [-200, 200];
                latlim = [-100, 100];
        end
    end

    function [xx,yy] = region_xxyy(region)
        [lonlim, latlim] = define_regions(region);
        [glon, glat] = geos_chem_centers;
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
        % 'seas' (southeast asia).
        
        % Parse the optional arguments
        fig = [];
        column_bool = false;
        p_max = [];
        p_mat = [];
        region = '';
        tp_mask = [];
        for a=1:numel(varargin)
            if isnumeric(varargin{a}) && isscalar(varargin{a})
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
            end
        end
        
        [lonlim, latlim] = define_regions(region);
        
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
            elseif isempty(tp_mask)
                E.badinput('A troposphere mask must be given if plotting a concentration')
            elseif isstruct(tp_mask) && ~strcmpi(tp_mask.fullName,'')
                E.badinput('The tropopause mask is a structure but does not seem to refer to the tropopause level')
            end
            
            % Set stratospheric boxes to be NaN
            if isstruct(tp_mask)
                tp_mask = trop_mask(size(dataStruct.dataBlock), tp_mask);
            end
            dataStruct.dataBlock(~tp_mask) = nan;
            
            % Average the concentrations over the requested levels. Set
            % boxes below the max pressure to be NaN.
            if ndims(dataStruct.dataBlock) ~= ndims(p_mat) || any(size(dataStruct.dataBlock) ~= size(p_mat))
                E.badinput('Pressure and concentration matrices must be the same size')
            end
            
            pp = p_mat <= p_max;
            dataStruct.dataBlock(~pp) = nan;
            data = squeeze(nanmean(dataStruct.dataBlock,3));
        end
        
        % criteria must be a logical matrix with the same lat/lon/time
        % dimensions as dataBlock or Columns in the input structure.
        if any(size(data) ~= size(criteria))
            E.badinput('Size of input data and criteria do not match')
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
        
        avg_data = nanmean(nanmean(data,3),2);
        
        
        
        % Create or select the figure as necessary
        if isempty(fig)
            figure;
        else
            figure(fig);
            hold on
        end
        
        plot(glon(xx)', avg_data);
        
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



end

