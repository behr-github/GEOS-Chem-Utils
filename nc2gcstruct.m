function [ GC ] = nc2gcstruct( ncfile )
%NC2GCSTRUCT Read a netCDF version of extracted GEOS-Chem output into a structure
%   For the lightning NOx lifetime paper, we used certain utility functions
%   that read GEOS-Chem binary punch (bpch) files into Matlab structures
%   with a certain format. To make that data available with the publication
%   in an open-source format, we saved the data we used as netCDF4 files.
%   This function will read those files and return a similar structure to
%   the one we used.

ni = ncinfo(ncfile);

% Find the variables that represent data, rather than coordinates, by
% assuming that coordinate variables are 1-D and data variables are 3- or
% 4-D
data_vars = false(1, numel(ni.Variables));
for a=1:numel(data_vars)
    data_vars(a) = numel(ni.Variables(a).Dimensions) > 1;
end

varnames = {ni.Variables.Name};
data_varnames = varnames(data_vars);

GC = repmat(struct('dataBlock',[],'dataUnit','','fullName','','tVec',[],'modelName','','modelRes',[],'lon',[],'lat',[]), 1, sum(data_vars));

for a=1:numel(GC)
    % Read in the common values and attributes first
    GC(a).modelName = ncreadatt(ni.Filename, '/', 'modelName');
    GC(a).modelRes = ncreadatt(ni.Filename, '/', 'modelRes');
    GC(a).tVec = ncread(ni.Filename, 'time');
    GC(a).lon = ncread(ni.Filename, 'lon');
    GC(a).lat = ncread(ni.Filename, 'lat');
    
    % Now the variable-specific values
    GC(a).dataBlock = ncread(ni.Filename, data_varnames{a});
    GC(a).dataUnit = ncreadatt(ni.Filename, data_varnames{a}, 'Units');
    % The variable names have been sanitized a bit from how they are read
    % directly from the .bpch files, in particular, excess dashes (-) and
    % dollar signs ($) were removed. Some functions that rely on matching
    % the variable name may need to be modified to take that into account.
    GC(a).fullName = data_varnames{a};
end


end

