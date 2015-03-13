function [ ctmData,mdlData,bpchInfo ] = readAllBPCHData( inputFile,tracerInfo,diagInfo,verbose,bruteForce,onlyInfo )
%READALLBPCHDATA Retrieves information from BPCH file
%   Extracts all data from a BPCH file as generated by GEOS-Chem. Files
%   which can be read are:
%       ctm.bpch:           Primary GC output file
%       tsYYYYMMDD.bpch:    ND49 timeseries output
%       restart.bpch:       GC restart files
%   Input arguments (default values in brackets, where appropriate):
%       inputFile:          Filepath of the BPCH file
%       tracerInfo:         Filepath of the tracerinfo.dat file
%                           (<simulation directory>/tracerinfo.dat)
%       diagInfo:           Filepath of the diaginfo.dat file
%                           (<simulation directory>/diaginfo.dat)
%       verbose:           	(T)/F: Show all warnings
%       bruteForce:         T/(F): Ignore non-fatal errors
%       onlyInfo:           T/(F): Retrieve only bpchInfo
%   Default tracerinfo and/or diaginfo files are used if:
%       i)   Relevant argument is given as 'default', or
%       ii)  Relevant argument is an empty array, or
%       iii) Relevant argument is not given
%   Model data is stored under ctmData.modelData, with these fields:
%       hRes:       2x1 model resolution array
%       name:       Cell containing the model name
%       hPolar:     Boolean. True if model has half-sized boxes at poles
%       c180:       Boolean. True if the first longitude is centred on 180
%       longC:      Longitude centres
%       latC:       Latitude centres
%       longE:      Longitude edges
%       latE:       Latitude edges
%       numAlts:    Number of altitude levels
%   Tracer data is stored under ctmData.{category ID}.{tracer ID} as follows:
%       name        Full tracer name
%       data:       NIxNJxNL array of data
%       dims:       [NI NJ NL] for the tracer
%       unit:       Tracer unit, e.g. 'ppbv'
%       start:      [I0 J0 L0] for the tracer
%       weight:     Tracer molecular weight
%       carbon:     Moles of carbon/mole of tracer
%       datenum:    Vector of date numbers corresponding to the simulation
%                   date of data storage
%   Category names are converted to 'machine-safe' strings, with 'C_'
%   placed in front. Therefore 'IJ-AVG-$' becomes 'C_IJ_AVG'. The same
%   operation is performed on tracer IDs to create usable structure fields.
%
%   Example:
%       ctmData.C_IJ_AVG.T_NOx.dims = [47 46 72]
%
%   Note: If only bpchInfo is requested, 

if nargout > 1
    getMdlData = true;
    acquiredPSurf = false;
else
    getMdlData = false;
end

if ~exist('onlyInfo','var')
    onlyInfo = false;
end
if ~exist('verbose','var')
    verbose = true;
end
if ~exist('bruteForce','var')
    bruteForce = false;
end

defaultTracerInfo = true;
if exist('tracerInfo','var') && ~isempty(tracerInfo) && ~strcmpi(tracerInfo,'default')
    defaultTracerInfo = false;
end

defaultDiagInfo = true;
if exist('diagInfo','var') && ~isempty(diagInfo) && ~strcmpi(diagInfo,'default')
    defaultDiagInfo = false;
end

if defaultTracerInfo
    % Use tracerinfo.dat from same folder as input file
    tracerInfo = inputFile;
    % Change all backslashes into slashes
    tracerInfo(regexp(tracerInfo,'\')) = '/';
    tracerEnd=regexp(tracerInfo,'/');
    tracerInfo=sprintf('%stracerinfo.dat',tracerInfo(1:tracerEnd(end)));
end

if defaultDiagInfo
    % Use diaginfo.dat from same folder as input file
    diagInfo = inputFile;
    % Change all backslashes into slashes
    diagInfo(regexp(diagInfo,'\')) = '/';
    diagEnd=regexp(diagInfo,'/');
    diagInfo=sprintf('%sdiaginfo.dat',diagInfo(1:diagEnd(end)));
end
    
% Check that files exist
if ~exist(inputFile,'file')
    error('BPCHRead:MissingTarget','BPCH input file not found.');
elseif ~exist(tracerInfo,'file')
    error('BPCHRead:MissingTracerInfo','Tracer data file not found.');
elseif ~exist(diagInfo,'file')
    error('BPCHRead:MissingDiagInfo','Diagnostics data file not found.');
end

% Open file
fileID = fopen(inputFile,'r','ieee-be');

% Get header information
%[fileType,~,rOK_A] = readFORTRANRecord(fileID,'*char',1);
[~,~,rOK_A] = readFORTRANRecord(fileID,'*char',1);
%fileType = strtrim(fileType');
[titleLine,~,rOK_B] = readFORTRANRecord(fileID,'*char',1);
titleLine = strtrim(titleLine');
testRead(rOK_A & rOK_B);

% Check that we recognise the title line
PSCDiagTitle = 'GEOS-CHEM Checkpoint File: Instantaneous PSC state (unitless)';
CSPECTitle = 'geos-chem checkpoint file: instantaneous species concentrations (#/cm3)';
knownTitlesLC={'geos-chem diag49 instantaneous timeseries';...      % Timeseries data
    'geos-chem binary punch file v. 2.0';...                        % Standard CTM output
	'geos-chem adj file: instantaneous adjoint concentrations'};    % Adjoint output

% Special case
% CSPEC/PSC diagnostic file info not yet sent to tracerinfo/diaginfo - handle manually
isPSCDiag = strcmpi(titleLine,PSCDiagTitle);
isCSPEC = strcmpi(titleLine,CSPECTitle);
if ~(isPSCDiag || isCSPEC || any(strcmpi(titleLine,knownTitlesLC)))
    warning('BPCHRead:UnknownFileType','File type ''%s'' not recognised; attempting to parse.',titleLine);   
end

% Store this position - will need to rewind to it
dataStartPos = ftell(fileID);

%% Read in data to establish main data structure

% First line of the header of the first datablock
[rOK,modelName,modelRes,halfPolar,centre180]=...
    readFixedFORTRANRecord(fileID,'*char',20,'*float32',2,'*int32',1,'*int32',1);

modelName = strtrim(modelName');

testRead(rOK);

longRes = modelRes(1);
latRes = modelRes(2);

% Horizontal resolution options
%   GCAP 4x5
%   GMAO 4x5
%   GMAO 2x2.5
%   GMAO 1x1.25
%   GMAO 1x1
%   GMAO 0.5x0.667
%   SE Asia Nested Grid
%   N America Nested Grid
%   Europe Nested Grid: 
%   Generic 1x1: 360x180 cells

% Determine the following:
%   latEVec:  Vector of latitude edges
%   latCVec:  Vector of latitude centres
%   longEVec: Vector of longitude edges
%   longCVec: Vector of longitude centres
%   numAlts:  Number of altitude levels

% Handle nested grid sims elsewhere
switch lower(modelName(1:4))
    case 'geos'
        if length(modelName) > 5
            % Reduced level model - only in the reduced level cases are the
            % modelnames more expressive than 'GEOS5' or 'GEOS3'
            reducedLev = true;
        else
            reducedLev = false;
        end
        switch lower(modelName(5))
            case '3'
                if reducedLev
                    numAlts = 48;
                else
                    numAlts = 30;
                end
            case '4'
                if reducedLev
                    numAlts = 30;
                else
                    numAlts = 55;
                end
            case '5'
                if reducedLev
                    numAlts = 47;
                else
                    numAlts = 72;
                end
        end
        % Determine horizontal cell numbers
        switch round(longRes*3)
            case 15
                % GMAO 4x5 grid
                longEVec = -182.5:5:177.5;
                longCVec = -180:5:175;
                latEVec  = [-90,-88:4:88,90];
                latCVec  = [-89,-86:4:86,89];
            case 8
                % GMAO 2x2.5 grid
                longEVec = -181.25:2.5:178.75;
                longCVec = -180:2.5:177.5;
                latEVec  = [-90,-89:2:89,90];
                latCVec  = [-89.5,-88:2:88,89.5];
            case 4
                longEVec = -180.625:1.25:179.375;
                longCVec = -180:1.25:178.75;
                latEVec  = [-90,-89.5:1:89.5,90];
                latCVec  = [-89.75,-89:1:89,89.75];
            case 3
                longEVec = -180.5:1:179.5;
                longCVec = -180:1:179;
                latEVec  = [-90,-89.5:1:89.5,90];
                latCVec  = [-89.75,-89:1:89,89.75];
            case 2
                oneThrd = (1/3);
                twoThrd = (2/3);
                longEVec = (-180-oneThrd):twoThrd:(179+twoThrd);
                longCVec = -180:twoThrd:(179+oneThrd);
                latEVec  = [-90,-89.75:0.5:89.75,90];
                latCVec  = [-89.875,-89.5:0.5:89,89.875];
        end
    case 'gcap'
        % GCAP 4x5
        longEVec = -182.5:5:177.5;
        longCVec = -180:5:175;
        latEVec  = -90:4:90;
        latCVec  = -88:4:88;
        numAlts = 23;
    otherwise
        error('readAllBPCHData:badRes','Grid ''%s'' not recognized',modelName);
end

%fprintf('MN: %s\n',modelName);

if isPSCDiag
    tID = {'STATE_PSC'};
    tName = {'PSC state'};
    tWeight = 1;
    tCarbon = 1;
    tNum = 1;
    tScale = 1;
    tUnit = {'-'};
    dOffset = 0;
    dName = {'IJ-PSC-$'};
elseif isCSPEC
    % Need to scan through CSPEC file to determine number of gases
    tID = {'STATE_PSC'};
    tName = {'PSC state'};
    tWeight = 1;
    tCarbon = 1;
    tNum = 1;
    tScale = 1;
    tUnit = {'molec/cm3/box'};
    dOffset = 0;
    dName = {'IJ-CHK-$'};
else
    % Read tracer database
    [tID,tName,tWeight,tCarbon,tNum,tScale,tUnit] = readTracerData(tracerInfo);

    % Read diagnostics information
    [dOffset,dName] = readDiagInfo(diagInfo);
end
numTracers = length(tWeight);
numDiags = length(dOffset);

%% Sanitise field names
% Some diagnostics and tracer have names that won't work as MATLAB fields
% Convert these to safe strings for field names
dNameSafe = dName;
for iDiag = 1:numDiags
    safeStr = char(dName{iDiag});
    if ~(isempty(regexp(safeStr,'-\$', 'once')) && isempty(regexp(safeStr,'=\$', 'once')))
        % Last two characters are -$ or =$
        safeStr = safeStr(1:end-2);
    end
    safeStr(regexp(safeStr,'\$')) = '';
    safeStr(regexp(safeStr,')')) = '';
    safeStr(regexp(safeStr,'(')) = '';
    safeStr(regexp(safeStr,'-')) = '_';
    safeStr(regexp(safeStr,' ')) = '_';
    safeStr(regexp(safeStr,'/')) = '_';
    safeStr(regexp(safeStr,'=')) = '_';
    % Ensure that the start of the category name is valid
    safeStr = ['C_' safeStr]; %#ok<AGROW>
    if ~isvarname(safeStr)
        if ~bruteForce
            error('BPCHRead:FieldSanity','Could not produce safe field name for diagnostic category %s (attempted ''%s'').',char(dName{iDiag}),safeStr);
        elseif verbose
            warning('BPCHRead:FieldSanity','Could not produce safe field name for diagnostic category %s (attempted ''%s''). Ignoring error.',char(dName{iDiag}),safeStr);
        end
    end
    dNameSafe(iDiag) = {safeStr};
end
tIDSafe = tID;
for iTracer = 1:numTracers
    safeStr = char(tID{iTracer});
    safeStr(regexp(safeStr,'\$')) = '';
    safeStr(regexp(safeStr,')')) = '';
    safeStr(regexp(safeStr,'(')) = '';
    safeStr(regexp(safeStr,'-')) = '_';
    safeStr(regexp(safeStr,' ')) = '_';
    safeStr(regexp(safeStr,'/')) = '_';
    safeStr(regexp(safeStr,'=')) = '_';
    % Ensure that the start of the category name is valid
    safeStr = ['T_' safeStr]; %#ok<AGROW>
    if ~isvarname(safeStr)
        if ~bruteForce
            error('BPCHRead:FieldSanity','Could not produce safe field name for tracer %s (attempted ''%s'').',char(tID{iTracer}),safeStr);
        elseif verbose
            warning('BPCHRead:FieldSanity','Could not produce safe field name for tracer %s (attempted ''%s''). Ignoring error.',char(tID{iTracer}),safeStr);
        end
    end
    tIDSafe(iTracer) = {safeStr};
end

bpchInfo.tracerIDs = tIDSafe;
bpchInfo.catIDs = dNameSafe;
%% Return if complete
if onlyInfo
    ctmData = 0;
    mdlData = 0;
    return
end

%% Establish main data structure

modelData = struct('hRes',[latRes longRes],'name',modelName,...
    'hPolar',halfPolar,'c180',centre180,'longE',longEVec,'latE',latEVec,...
    'longC',longCVec,'latC',latCVec,'numAlts',numAlts);
tracerStruct = struct('data',0,'dims',[],...
    'start',[],'weight',0,'carbon',0,'unit','','datenum',0);
ctmData.modelData = modelData;

%% Read datablocks
% Rewind to start of data (returns 0 if successful, -1 otherwise)
if fseek(fileID,dataStartPos,-1)
    error('Could not rewind to start of data blocks.');
end
fileComplete = false;

% Run a first pass - don't retrieve any data, just scan through the file to
% determine its structure
% Maximum possible number of categories determined by reading diagInfo
catStored = false(numDiags,numTracers);
catLength = zeros(numDiags,1);
currLoc = ftell(fileID);
fseek(fileID,0,1);
endLoc = ftell(fileID);
fseek(fileID,currLoc,-1);
catMissing = false;
while ~fileComplete
    % Header lines - formatting taken direct from website
    %[ROk1,MName,MRes,HPolar,C180]...
        %=readFixedFORTRANRecord(fileID,'*char',20,'*float32',2,'*int32',1,'*int32',1);
    % First line: model name, model resolution, and data concerning
    % formatting. We are assuming that this data will be the same
    % throughout the file, and has already been acquired.
    readFORTRANRecord(fileID,'seekpast');
    % Second line: information about the tracer itself.
    %   Diagnostic category name
    %   Tracer number
    %   Unit string (e.g. ppbv)
    %   Starting date (tau value - hours since 1/1/1985)
    %   Ending date
    %   Unused string
    %   Dimensions - NI, NJ, NL, I0, J0, L0
    %   Length of data block in bytes
    [readOK,diagCat,tracer,~,~,~,~,dataDim,~]...
        =readFixedFORTRANRecord(fileID,'*char',40,'*int32',1,'*char',40,...
        '*float64',1,'*float64',1,'*char',40,'*int32',6,'*int32',1);
    testRead(readOK);
    diagCat = strtrim(diagCat');
    catIndex = find(strcmp(diagCat,dName));
    catName = char(dNameSafe(catIndex));
    tracer = tracer + dOffset(catIndex);
    tracIndex = find(tNum==tracer);
    tracName = char(tIDSafe(tracIndex));
    numDims = max(sum(dataDim~=1),1);
    if isempty(catName)
        % No categories
        if ~catMissing
            if verbose
                warning('BPCHRead:NoCategories','No categories found; storing all tracers under category ''data''.');
            end
            % Reduce size of category matrices to 1 category
            catStored = catStored(1,:);
            catLength = 0;
        end
        catMissing = true;
        catName = 'data';
        catIndex = 1;
    end
    if isempty(tracName)
        % Missing tracer in tracerinfo.dat
        if verbose && bruteForce
            warning('BPCHRead:MissingTracer','Tracer %i not recorded in tracerinfo.dat - skipping data.',tracer);
        else
            error('BPCHRead:MissingTracer','Tracer %i not recorded in tracerinfo.dat - aborting.',tracer);
        end
    else
        if ~catStored(catIndex,tracIndex)
            % First instance of this category/tracer combination
            catStored(catIndex,tracIndex) = true;
            tempStruct = tracerStruct;
            tempStruct.name = tName{tracIndex};
            tempStruct.dims = dataDim(1:numDims)';
            tempStruct.unit = tUnit{tracIndex};
            tempStruct.weight = tWeight(tracIndex);
            tempStruct.carbon = tCarbon(tracIndex);
            ctmData.(catName).(tracName) = tempStruct;
        end
        % Will need to divide the category length by the number of tracers
        % found for that category
        catLength(catIndex) = catLength(catIndex) + 1;
    end
    % Skip actual data block
    readFORTRANRecord(fileID,'seekpast');
    %tempBlock   =reshape(readFORTRANRecord(fileID,'*float32',4),dataDim(1:3)');
    currLoc = ftell(fileID);
    locDelta = endLoc - currLoc;
    if locDelta == 0
        fileComplete = true;
    end
end

% Divide stored category lengths by the number of tracers found
catLength = catLength./sum(catStored,2);
catLength(isnan(catLength)) = 0;
catLength(isinf(catLength)) = 0;

% Now read in real data
if fseek(fileID,dataStartPos,-1)
    error('Could not rewind to start of data blocks.');
end
fileComplete = false;
catComplete = ~catLength;
tracAcquired = zeros(size(catStored));

% Rewind to beginning of datablocks
fseek(fileID,dataStartPos,-1);
while ~fileComplete
    % Header lines - formatting taken direct from website
    %[ROk1,MName,MRes,HPolar,C180]...
        %=readFixedFORTRANRecord(fileID,'*char',20,'*float32',2,'*int32',1,'*int32',1);
    % First line: model name, model resolution, and data concerning
    % formatting. We are assuming that this data will be the same
    % throughout the file, and has already been acquired.
    readFORTRANRecord(fileID,'seekpast');
    % Second line: information about the tracer itself.
    %   Diagnostic category name
    %   Tracer number
    %   Unit string (e.g. ppbv)
    %   Starting date (tau value - hours since 1/1/1985)
    %   Ending date
    %   Unused string
    %   Dimensions - NI, NJ, NL, I0, J0, L0
    %   Length of data block in bytes
    [readOK,diagCat,tracer,~,tauStart,tauEnd,~,dataDim,~]...
        =readFixedFORTRANRecord(fileID,'*char',40,'*int32',1,'*char',40,...
        '*float64',1,'*float64',1,'*char',40,'*int32',6,'*int32',1);
    testRead(readOK);
    diagCat = strtrim(diagCat');
    if ~catMissing
        catIndex = find(strcmp(diagCat,dName));
        catName = char(dNameSafe(catIndex));
    else
        catIndex = 1;
        catName = 'data';
    end
    tracer = tracer + dOffset(catIndex);
    tracIndex = find(tNum==tracer);
    tracName = char(tIDSafe(tracIndex));
    if isempty(tracName)
        % Missing tracer in tracerinfo.dat
        readFORTRANRecord(fileID,'seekpast');
    else
        timeIndex = tracAcquired(catIndex,tracIndex) + 1;
        tracAcquired(catIndex,tracIndex) = timeIndex;
        numDims = max(sum(dataDim~=1),1);
        if timeIndex == 1
            % First instance of this category
            ctmData.(catName).(tracName).datenum = zeros(catLength(catIndex),2);
            ctmData.(catName).(tracName).data = zeros([dataDim(1:numDims)',catLength(catIndex)]);
        end
        ctmData.(catName).(tracName).datenum(timeIndex,:) = tauToDate([tauStart tauEnd]);
        indexArray = cell(1,numDims+1);
        indexArray(1:numDims) = {':'};
        indexArray(end) = {timeIndex};
        %{
        tempData = readFORTRANRecord(fileID,'*float32',4);
        tempData = reshape(tempData,dataDim(1:numDims)');
        tempData = tempData .* tScale(tracIndex);
        ctmData.(catName).(tracName).data(indexArray{:}) = tempData;
        %}
        if getMdlData && ~acquiredPSurf && strcmpi('T_PSURF',tracName)
            % Pressure surface data found
            acquiredPSurf = true;
            PSCat = catName;
        end
        ctmData.(catName).(tracName).data(indexArray{:}) = tScale(tracIndex) .* reshape(readFORTRANRecord(fileID,'*float32',4),dataDim(1:numDims)');
    end
    currLoc = ftell(fileID);
    locDelta = endLoc - currLoc;
    if ((locDelta == 0) || ~(sum(sum(~catComplete))))
        fileComplete = true;
    end
end

if isPSCDiag
    % Correct to remove 0.1 offset
    ctmData.C_IJ_PSC.T_STATE_PSC.data = floor(ctmData.C_IJ_PSC.T_STATE_PSC.data);
end

% Create model data structure
if getMdlData
    mdlData.lonVals = longEVec;
    mdlData.latVals = latEVec;
    mdlData.zVals = [];
    mdlData.zUnit = 'hPa';
    gridArea = GEOSGridArea(longEVec,latEVec);
    mdlData.gridArea = gridArea;
    % Pressure levels problematic - first check if the full surface data
    % array is available
    if acquiredPSurf
        mdlData.pArray = ctmData.(PSCat).T_PSURF.data;
        mdlData.zUnit = ctmData.(PSCat).T_PSURF.unit;
        meanPLevs = (mean(mean(mdlData.pArray,4),1));
        % Pressure levels are now averaged by longitude and time, but
        % averaging by latitude is more difficult
        gridArea = repmat(gridArea(1,:),[1,1,size(meanPLevs,3),1]);
        meanPLevs = (meanPLevs.*gridArea)./sum(gridArea(1,:,1,1));
        meanPLevs = squeeze(sum(meanPLevs,2));
        % We know that the top level is always 0.01hPa
        if meanPLevs(end) > 0.01001
            meanPLevs = [meanPLevs(:);0.01];
        end
        mdlData.zVals = meanPLevs;
    end
end

% ----- Future work: check for multiple entries with the same timestamp,
% category and tracer -----

fclose(fileID);

    function []=testRead(readOK)
        if ~readOK
            fclose(fileID);
            error('Read error in file %s.\n',inputFile);
        end
    end

end

function MATDate = tauToDate(tauDate)
% 725008 = datenum(1985,1,1,0,0,0)
    MATDate = 725008 + (tauDate./24);
end

