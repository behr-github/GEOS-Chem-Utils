function [  ] = gcstruct2ncdf( ncfile, varargin )
%GCSTRUCT2NCFILE Save a GEOS-Chem structure as a netCDF file.
%   Detailed explanation goes here

[glon,glat] = geos_chem_centers('2x25');

xx = strcmpi(varargin,'overwrite');
if any(xx)
    overwrite = true;
    varargin(xx) = [];
else
    overwrite = false;
end

if exist(ncfile, 'file')
    if overwrite || ask_yn(sprintf('File %s exists. Overwrite?', ncfile))
        delete(ncfile);
    else
        error('io:file_exists', 'File %s exists, aborting', ncfile)
    end
end          

% Fields that should be the same in all input structures
common_fields = {'modelName', 'tEdge', 'tVec', 'modelRes'};
first_struct = true;
for a=1:numel(varargin)
    for b=1:numel(varargin{a})
        GC = varargin{a}(b);
        if first_struct
            % Write the common fields, as global attributes or variables as
            % needed
            write_var_with_atts(ncfile,'lon',glon,make_dims(glon, GC.tVec,GC.tEdge), 'degrees', 'Longitude of GEOS-Chem grid cell center (west is negative)');
            write_var_with_atts(ncfile,'lat',glat,make_dims(glat, GC.tVec,GC.tEdge), 'degrees', 'Latitude of GEOS-Chem grid cell center (south is negative)');
            write_var_with_atts(ncfile,'time',GC.tVec - datenum('1985-01-01'),make_dims(GC.tVec,GC.tVec,GC.tEdge), 'days since midnight 1 Jan 1985', 'Model date as a number of days since midnight 1 Jan 1985');
            write_var_with_atts(ncfile,'time_edge',GC.tEdge - datenum('1985-01-01'),make_dims(GC.tEdge,GC.tVec,GC.tEdge),'days since midnight 1 Jan 1985','Edge of model averaging periods (used only for variables with adjacent averaging periods)');
            ncwriteatt(ncfile,'/','modelName',GC.modelName);
            ncwriteatt(ncfile,'/','modelRes',GC.modelRes);
            FirstGC = GC;
            first_struct = false;
        else
            for c=1:numel(common_fields)
                if ~isequaln(GC.(common_fields{c}), FirstGC.(common_fields{c}))
                    error('common_field:not_equal', 'Common field %s not equal in structure %d index %d and first structure', common_fields{b}, a, b);
                end
            end
        end
        gcvar = regexprep(GC.fullName, '\W', ''); % Remove everything but a-z, 0-9, and _
        dims_cell = make_dims(GC.dataBlock, GC.tVec, GC.tEdge);
        write_var_with_atts(ncfile, gcvar, GC.dataBlock, dims_cell, GC.dataUnit, GC.description);
    end
end

    function dims = make_dims(dataBlock, tvec, tedge)
        dims = cell(1, 2*ndims(dataBlock));
        sz = size(dataBlock);
        sz_matches = {numel(glon), 'lon';
                      numel(glat), 'lat';
                      numel(tvec), 'time';
                      numel(tedge), 'time_edge';
                      47, 'model_level'};
        dim_len_1 = false(size(dims));
        for i=1:numel(sz)
            i_name = (i-1)*2 + 1;
            i_length = i_name + 1;
            if sz(i) == 1
                dim_len_1([i_name, i_length]) = true;
                continue
            end
            xx = sz(i) == cat(1,sz_matches{:,1});
            if sum(xx) == 1
                dims{i_name} = sz_matches{xx,2};
                dims{i_length} = sz_matches{xx,1};
            else
                error('ncdim:cannot_id','Cannot identify dimension length of %d', sz(i));
            end
        end
        
        dims = dims(~dim_len_1);
    end

end

function write_var_with_atts(ncfile, gcvar, val, dims, units, description)
if ~isnumeric(val)
    error('ncval:not_numeric', 'VAL must be a numeric type');
end
val = single(val);
nccreate(ncfile, gcvar, 'Dimensions', dims, 'Datatype', 'single');
ncwrite(ncfile, gcvar, val);
ncwriteatt(ncfile, gcvar, 'Units', units);
ncwriteatt(ncfile, gcvar, 'Description', description);
end