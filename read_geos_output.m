function [ Data ] = read_geos_output(  )
%read_geos_output Uses the functions written by Sebastion Eastham to pull
%in variables from geos chem output (binary punch) files.
%  	The read BPCH function require a goodly number of inputs; this uses the
%  	MATLAB UI features to make it slightly more user friendly.  It will ask
%  	for the category (diagnostic) and tracer name - these should be the
%  	"Matlab sanitized" versions.
%
%   This function assumes that the tracerinfo.dat and diaginfo.dat files
%   are in the same folder as the output file of interest. If that is not
%   true, an error will occur.

% Initialize the error handling class
E = JLLErrors;

% Add the folder with the BPCH functions. In should be in the folder with
% this function.
addpath('~/Documents/MATLAB/GEOS_Chem_Utils/BPCH_Functions');

% Get the output file from the user, and use the directory to find the
% tracerinfo and diaginfo files.

getfile_title = 'Select the output file to read. tracerinfo and diaginfo should be in the same folder.';
[output_filename, pathname] = uigetfile('*',getfile_title);
if output_filename == 0;
    error(E.userCancel);
end

if ~exist(fullfile(pathname,'tracerinfo.dat'),'file') || ~exist(fullfile(pathname,'diaginfo.dat'),'file')
    error(E.callError('dat_files_missing','The tracerinfo.dat and diaginfo.dat files must in the same folder as the output file.'));
end

inputFile = fullfile(pathname,output_filename);
tracerFile = fullfile(pathname,'tracerinfo.dat');
diagFile = fullfile(pathname,'diaginfo.dat');

Data = struct('dataBlock', [], 'dataUnit', [], 'fullName', [], 'fullCat', [], 'tVec', [], 'modelName', [], 'modelRes', [], 'dataScale', [], 'molMass', [], 'tEdge', []);
input_title = 'Enter the %s, or cancel to exit and return.';
D=0;
% Ask the user for the category and tracer names.  Loop until the user
% cancels or enters an empty string for category or tracer; this allows the
% user to enter multiple variables.  Note that inputdlg returns an empty
% cell array, i.e. {}, if the user cancels and a cell array with an empty
% string, i.e. {''}, if the user hits OK without entering anything.  We'll
% use this to check if the user wants to append the data useful for dealing
% with satellites, such as pressure, boxheight, etc.
first_time = true;
while true
    user_input = inputdlg(sprintf(input_title,'category'));
    if isempty(user_input);
        append_sat = false;
        break
    elseif isempty(user_input{1})
        append_sat = true;
        break
    end
    
    category = user_input{1};
    
    user_input = inputdlg(sprintf(input_title,'tracer'));
    if isempty(user_input)
        append_sat = false;
        break
    elseif isempty(user_input{1})
        append_sat = true;
        break
    end
    
    tracer = user_input{1};
    
    % Change the message for successive variables.
    if first_time
        first_time = false;
        input_title = 'Enter another %s, enter a blank field to append the variables NAIR, PSURF,and BXHEIGHT, or cancel to return immediately.';
    end
    
    % Call the function to read the BPCH file
    try
        [ dataBlock, dataUnit, fullName, fullCat, tVec, modelName, modelRes, dataScale, molMass, tEdge ] = readBPCHSingle(inputFile,category,tracer,tracerFile,diagFile);
        % Increment the index of the return structure. Placing this after
        % the read function means that if there is an caught error in it,
        % we won't unnecessarily advance the counter.
        D=D+1;
        Data(D).dataBlock = dataBlock;
        Data(D).dataUnit = dataUnit;
        Data(D).fullName = fullName;
        Data(D).fullCat = fullCat;
        Data(D).tVec = tVec;
        Data(D).modelName = modelName;
        Data(D).modelRes = modelRes;
        Data(D).dataScale = dataScale;
        Data(D).molMass = molMass;
        Data(D).tEdge = tEdge;
    catch err
        % Catch errors resulting from the user entering a wrong category or
        % tracer.  This will present a message box that waits until the
        % user presses OK before continuing.
        if strcmp(err.identifier,'BPCHRead:UnknownTracer')
            tr_or_cat = 'tracer';
            couldnt_find = tracer;
        elseif strcmp(err.identifier,'BPCHRead:UnknownCategory')
            tr_or_cat = 'category';
            couldnt_find = category;
        else
            rethrow(err);
        end
        box_msg = sprintf('The %s %s was invalid. This category & tracer will not be read.',tr_or_cat,couldnt_find);
        box_title = sprintf('Could not find %s',tr_or_cat); 
        uiwait(msgbox(box_msg,box_title,'error','modal'));
    end
end

% Get the NAIR, PSURF, and BXHEIGHT variables, if the user hit "OK" on a
% blank input
if append_sat
    categories = {'BXHGHT','PEDGE','BXHGHT','TR_PAUSE'};
    tracers = {'NAIR','PSURF','BXHEIGHT','TP_LEVEL'};
    for a=1:numel(categories);
        fprintf('Now retrieving %s/%s\n',categories{a},tracers{a});
        D = D+1;
        [ dataBlock, dataUnit, fullName, fullCat, tVec, modelName, modelRes, dataScale, molMass, tEdge ] = readBPCHSingle(inputFile,categories{a},tracers{a},tracerFile,diagFile);
        Data(D).dataBlock = dataBlock;
        Data(D).dataUnit = dataUnit;
        Data(D).fullName = fullName;
        Data(D).fullCat = fullCat;
        Data(D).tVec = tVec;
        Data(D).modelName = modelName;
        Data(D).modelRes = modelRes;
        Data(D).dataScale = dataScale;
        Data(D).molMass = molMass;
        Data(D).tEdge = tEdge;
    end
end


end

