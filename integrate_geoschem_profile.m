function [ StructIn ] = integrate_geoschem_profile( StructIn, unit_conv, time_indicies )
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
%   There is also an optional argument:
%       
%       3) A description of which time indicies to use.  This can either be
%       the literal 4th dimension indicies (as a row vector) or a 2-element
%       cell array with start and end dates as strings that the datenum()
%       function will recognize.



% Declaring error handler class

E = JLLErrors;

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


% Check that if time indicies are given, they make up a or a cell array
if nargin > 2 && ~(isvector(time_indicies) || iscell(time_indicies))
    error(E.badinput('''time_indicies'' must be a row vector or 2-element cell array, if specified'))
elseif nargin > 2 && iscell(time_indicies) && numel(time_indicies) ~= 2
    error(E.badinput('''time_indicides'' must have two elements if passed as a cell array'));
elseif nargin > 2 && ~isrow(time_indicies) && isnumeric(time_indicies)
    error(E.badinput('''time_indicides'' must be a row vector'));
end

% Find the number density of air, pressure levels, and box height fields.
% If one of them is not present, throw an error.
nair_ind = 0;
psurf_ind = 0;
bxheight_ind = 0;
for i=1:numel(StructIn)
    if strcmp(StructIn(i).fullName,'Number density of air')
        nair_ind = i;
    elseif strcmp(StructIn(i).fullName,'Surface pressure')
        psurf_ind = i;
    elseif strcmp(StructIn(i).fullName,'Grid box height')
        bxheight_ind = i;
    elseif strcmp(StructIn(i).fullName,'Tropopause level')
        tplevel_ind = i;
    end
end

if nair_ind == 0
    error(E.badinput('The input structure must contain a fullName field with the value ''Number density of air'''));
elseif psurf_ind == 0
    error(E.badinput('The input structure must contain a fullName field with the value ''Surface pressure'''));
elseif bxheight_ind == 0
    error(E.badinput('The input structure must contain a fullName field with the value ''Grid box height'''));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VARIABLE PREP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy the nair, pressure, and box height values
nair = StructIn(nair_ind).dataBlock;
psurf = StructIn(psurf_ind).dataBlock;
bxheight = StructIn(bxheight_ind).dataBlock .* 100; % Box height is given in meters, we need cm
tplevel = floor(StructIn(tplevel_ind).dataBlock); % GEOS-Chem gives partial levels for its tropopause, we want only levels entirely in the troposphere
s = size(nair);

% If time indicies is a vector of indicies, do nothing. If it wasn't
% passed, assume that all times should be integrated.  If a cell array was
% given, find all times between the start and end dates.

if ~exist('time_indicies','var')
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
    end
    
    % Prepare a matrix with dimensions lon x lat x time for the integrated
    % columns.
    columns = zeros(s(1),s(2),s(4));
    
    % Iterate through the timesteps, even if there's just one
    for a=time_indicies;

        nair_subset = nair(:,:,:,a);
        p_edges_subset = psurf(:,:,:,a);
        bxheight_subset = bxheight(:,:,:,a);
        species_subset = StructIn(b).dataBlock(:,:,:,a);
        tplevel_subset = tplevel(:,:,a);
        
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
                
                % Restrict the data vectors to those indicies wholly in the
                % troposphere
                bxheight_vec = bxheight_vec(1:tplevel_val);
                species_vec = species_vec(1:tplevel_val);
                
                % Integrate the concentration using the trapzoid rule
                %columns(i,j,a) = trapz(bxheight_vec, species_vec);
                columns(i,j,a) = sum(bxheight_vec .* species_vec);
            end
        end
    end
    StructIn(b).Columns = columns;
end

end

