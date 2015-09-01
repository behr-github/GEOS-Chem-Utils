function [ StructIn ] = integrate_geoschem_profile( StructIn, unit_conv, time_indicies, cloud_column_sel )
%integrate_geoschem_profile Integrates GEOS-Chem output tracers into columns
%   While GEOS-Chem does have the capability of outputting instantaneous
%   columns for comparison to satellite overpasses, if you didn't output
%   that diagnostic or want a monthly averaged column, this function will
%   do that. Takes two required arguments:
%
%       1) A structure obtained by using read_geos_output. This should
%       automatically contain the "tracers" NAIR, BXHEIGHT, TP_LEVEL and
%       PSURF. 
%       2) The unit conversion from parts-per-(whatever) to
%       parts-per-part. for the mixing ratios, e.g. 1e-9 for ppb
%
%   There are also two optional argument:
%       
%       3) A description of which time indicies to use.  This can either be
%       the literal 4th dimension indicies (as a row vector) or a 2-element
%       cell array with start and end dates as strings that the datenum()
%       function will recognize. Pass as an empty matrix to use all time
%       indicies and enter the fourth variable.
%
%       4) An integer (defaults to 0) whether to only integrate above
%       clouds. 1 will just integrate each cell above the cloud top. 2 will
%       do an average of clear and cloudy column weighted by the cloud
%       fraction.



% Declaring error handler class

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING AND CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for 2 to 4 inputs
narginchk(2,4);

% Check that the first argument is a structure
if ~isstruct(StructIn)
    error(E.badinput('Input ''StructIn'' must be a structure'));
end

% Check that the unit conversion is a scalar
if ~isscalar(unit_conv)
    error(E.badinput('Input ''unit_conv'' must be a scalar'));
end

% Make sure the fourth argument is a scalar. Set it to 0 if not passed.
if nargin < 4
    cloud_column_sel = 0;
end
if ~isscalar(cloud_column_sel) || cloud_column_sel < 0 || cloud_column_sel > 2
    E.badinput('The only allowed options for cloud_column are 0 (ignore clouds) 1 (only above clouds) 2 (VCD weighted by cloud fraction)');
end

% Check that if time indicies are given, they make up a or a cell array
if nargin > 2 && ~(isvector(time_indicies) || iscell(time_indicies) || isempty(time_indicies))
    error(E.badinput('''time_indicies'' must be a row vector or 2-element cell array, if specified'))
elseif nargin > 2 && iscell(time_indicies) && numel(time_indicies) ~= 2
    error(E.badinput('''time_indicides'' must have two elements if passed as a cell array'));
elseif nargin > 2 && ~isrow(time_indicies) && isnumeric(time_indicies) && ~isempty(time_indicies)
    error(E.badinput('''time_indicides'' must be a row vector'));
end

% Find the number density of air, pressure levels, and box height fields.
% If one of them is not present, throw an error.
nair_ind = 0;
psurf_ind = 0;
bxheight_ind = 0;
for i=1:numel(StructIn)
    if ~isempty(regexp(StructIn(i).fullName,'Number density of air', 'once')) % use regexp here to let the field be "Number density of air" or "Number density of air (calculated)"
        nair_ind = i;
    elseif strcmp(StructIn(i).fullName,'Surface pressure') || ~isempty(regexp(StructIn(i).fullName,'PSURF', 'once'))
        psurf_ind = i;
    elseif strcmp(StructIn(i).fullName,'Grid box height') || ~isempty(regexp(StructIn(i).fullName,'BXHEIGHT', 'once'))
        bxheight_ind = i;
    elseif strcmp(StructIn(i).fullName,'Tropopause level')
        tplevel_ind = i;
    end
end

if nair_ind == 0
    error(E.badinput('The input structure must contain a fullName that includes ''Number density of air'''));
elseif psurf_ind == 0
    error(E.badinput('The input structure must contain a fullName field with the value ''Surface pressure'' or include the word ''PSURF'''));
elseif bxheight_ind == 0
    error(E.badinput('The input structure must contain a fullName field with the value ''Grid box height'' or include the word ''BXHEIGHT'''));
end

% Find the cloud pressure field and cloud fraction field - if they're
% needed and not found throw an error.
cld_top_ind = 0;
cld_frac_ind = 0;
for i=1:numel(StructIn)
    if strcmp(StructIn(i).fullName,'GMAO CLDTOP field')
        cld_top_ind = i;
    elseif strcmp(StructIn(i).fullName,'GMAO CLDFRC field')
        cld_frac_ind = i;
    end
end

if cld_top_ind == 0 && cloud_column_sel > 0
    E.badinput('The input structure must contain a fullName field ''GMAO CLDTOP field'' to use cloud_column > 0');
elseif cld_frac_ind == 0 && cloud_column_sel > 1
    E.badinput('The input structure must contain a fullName field ''GMAO CLDFRC field'' to use cloud_column > 1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VARIABLE PREP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy the nair, pressure, and box height values
nair = StructIn(nair_ind).dataBlock;
psurf = StructIn(psurf_ind).dataBlock;
bxheight = StructIn(bxheight_ind).dataBlock .* 100; % Box height is given in meters, we need cm
tplevel = floor(StructIn(tplevel_ind).dataBlock); % GEOS-Chem gives partial levels for its tropopause, we want only levels entirely in the troposphere

if cloud_column_sel > 0
    cldtop = floor(StructIn(cld_top_ind).dataBlock); % This will again be in levels (seriously GEOS-Chem?)
end
if cloud_column_sel > 1
    cldfrac = StructIn(cld_frac_ind).dataBlock;
end
s = size(nair);

% If time indicies is a vector of indicies, do nothing. If it wasn't
% passed, assume that all times should be integrated.  If a cell array was
% given, find all times between the start and end dates.

if ~exist('time_indicies','var') || isempty(time_indicies);
    time_indicies = 1:numel(StructIn(nair_ind).tVec);
elseif iscell(time_indicies)
    dates = datestr(time_indicies);
    xx = StructIn(nair_ind).tVec > min(dates) && StructIn(nair_ind).tVec < max(dates);
    time_indicies = find(xx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALCULATE COLUMNS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iterate through all the species in the structure, but skip the auxiliary
% data.
for b=1:numel(StructIn)
    if b==nair_ind; continue;
    elseif b==psurf_ind; continue;
    elseif b==bxheight_ind; continue;
    elseif b==tplevel_ind; continue
    elseif b==cld_top_ind; continue;
    elseif b==cld_frac_ind; continue;
    elseif isempty(regexp(StructIn(b).dataUnit, 'pp.', 'once')); continue; % skips fields that don't have a mixing ratio unit
    end
    
    % Prepare a matrix with dimensions lon x lat x time for the integrated
    % columns.
    columns = zeros(s(1),s(2),numel(time_indicies));
    
    % Iterate through the timesteps, even if there's just one
    for a=time_indicies;
        if DEBUG_LEVEL > 0; fprintf('Calculating for timestep %d\n',a); end
        nair_subset = nair(:,:,:,a);
        p_edges_subset = psurf(:,:,:,a);
        bxheight_subset = bxheight(:,:,:,a);
        species_subset = StructIn(b).dataBlock(:,:,:,a);
        tplevel_subset = tplevel(:,:,a);
        if cloud_column_sel > 0
            cldtop_subset = cldtop(:,:,a);
        end
        if cloud_column_sel > 1
            cldfrac_subset = cldfrac(:,:,a);
        end
        
        % Convert the species input from mixing ratio to number density
        %              = (mixing ratio) * (num density of air) * (m^3 / cm^3) * (parts-per-part / parts-per-{m,b,tr}illion)
        species_subset = species_subset .* nair_subset .* 1e-6 .* unit_conv;
        
        % Now go through each column and integrate
        
        for i=1:s(1)
            for j=1:s(2)
                % Reshape the altitude and concentration matrices into vectors
                %bxheight_vec = cumsum(squeeze(bxheight_subset(i,j,:))); % trapz needs x-coordinates, but bxheight is the actual height of each box
                bxheight_vec = squeeze(bxheight_subset(i,j,:));
                species_vec = squeeze(species_subset(i,j,:));
                p_vec = squeeze(p_edges_subset(i,j,1:end-1)); % There is always one more pressure edge than box
                tplevel_val = tplevel_subset(i,j);
                
                if cloud_column_sel > 0
                    cldtop_val = cldtop_subset(i,j);
                end
                if cloud_column_sel > 1
                    cldfrac_val = cldfrac_subset(i,j);
                end
                
                % Restrict the data vectors to those indicies wholly in the
                % troposphere
                bxheight_vec = bxheight_vec(1:tplevel_val);
                species_vec = species_vec(1:tplevel_val);
                
                % Integrate the concentration using the centerpoint rule
                %columns(i,j,a) = trapz(bxheight_vec, species_vec);
                clear_column = sum(bxheight_vec .* species_vec);
                
                % Now do the same from the cloud top, if we're supposed to 
                if cloud_column_sel > 0
                    bxheight_vec_cld = bxheight_vec(cldtop_val:end);
                    species_vec_cld = species_vec(cldtop_val:end);
                    cloud_column= sum(bxheight_vec_cld .* species_vec_cld);
                end
                
                % Now we will determine how to output the columns.
                % cloud_column_sel = 0 means give the whole column from
                % ground to tropopause. cloud_column_sel = 1 means just
                % give the column from cloud top to tropopause.
                % cloud_column_sel = 2 means weight the last two by the
                % cloud fraction - this is approximating the OMI retrieval
                % with no ghost columns.
                switch cloud_column_sel
                    case 0
                        columns(i,j,a) = clear_column;
                    case 1
                        columns(i,j,a) = cloud_column;
                    case 2
                        columns(i,j,a) = clear_column * (1 - cldfrac_val) + cloud_column * cldfrac_val;
                end
            end
        end
    end
    StructIn(b).Columns = columns;
end

end

