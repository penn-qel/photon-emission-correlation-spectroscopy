function TTTRdata = TTTR_import_PTU(filenames)
% Function to import TTTR data from PicoQuant .ptu files
%
% Inputs:
%   filenames: string or cell array of strings giving paths to .ptu files
% 
% Output: TTTRdata structure containig the fields:
%   "Name": Filename without the .ptu extension
%   "MagicNumber": Code that identifies type of PicoQuant File
%   "Version": Version of file
%   "Tags": Structure giving data stored in header (called "tags" by PicoQuant)  
%   "Events": Array of data encoding channel and time in .ptu format (uint32)
%
% If multiple files are requested, they are loaded into separate elements
% of the TTTRdata structure.
%
% Important note: This function simply reads the contents of the .ptu file(s) into
% memory. A warning is issued if the combined size of these files exceeds 1
% GB (or the available system memory).  If the files are very large and no
% detailed analysis is required it would be better to process them on the
% fly.
%
% This code is modified from the "Read_PTU.m" demo provided by
% PicoQuant
% 

%% some constants
    tyEmpty8      = hex2dec('FFFF0008');
    tyBool8       = hex2dec('00000008');
    tyInt8        = hex2dec('10000008');
    tyBitSet64    = hex2dec('11000008');
    tyColor8      = hex2dec('12000008');
    tyFloat8      = hex2dec('20000008');
    tyTDateTime   = hex2dec('21000008');
    tyFloat8Array = hex2dec('2001FFFF');
    tyAnsiString  = hex2dec('4001FFFF');
    tyWideString  = hex2dec('4002FFFF');
    tyBinaryBlob  = hex2dec('FFFFFFFF');
    % RecordTypes
%     rtPicoHarpT3     = hex2dec('00010303');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $03 (PicoHarp)
%     rtPicoHarpT2     = hex2dec('00010203');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $03 (PicoHarp)
%     rtHydraHarpT3    = hex2dec('00010304');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $04 (HydraHarp)
%     rtHydraHarpT2    = hex2dec('00010204');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $04 (HydraHarp)
%     rtHydraHarp2T3   = hex2dec('01010304');% (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
%     rtHydraHarp2T2   = hex2dec('01010204');% (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $04 (HydraHarp)
%     rtTimeHarp260NT3 = hex2dec('00010305');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $05 (TimeHarp260N)
%     rtTimeHarp260NT2 = hex2dec('00010205');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $05 (TimeHarp260N)
%     rtTimeHarp260PT3 = hex2dec('00010306');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $06 (TimeHarp260P)
%     rtTimeHarp260PT2 = hex2dec('00010206');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $06 (TimeHarp260P)


%% Check files and file sizes

if iscell(filenames)
    nfiles = length(filenames);
else
    nfiles = 1;
    filenames = {filenames}; % form into cell array
end

filesizes = NaN(1,nfiles);

for ff=1:nfiles
    fdata = dir(filenames{ff});
    filesizes(ff) = fdata.bytes;
end
%% Skipping this step cause i don't want to figure out how to adapt it for Mac.  Is that foolish? maybe. we shall see. -RF
% Check total file size against available memory

%[~,sys] = memory;
%datasize = sum(filesizes);

%if datasize>sys.PhysicalMemory.Available
%    error('Total size of .ptu files (%u MB) exceeds available system memory (%u MB)!',1e-6*datasize,1e-6*sys.PhysicalMemory.Available);
%end

%if datasize>1e9
%    choice = questdlg(sprintf('Total file size to import is %u MB.  Do you want to continue?',datasize*1e-6),...
%        'Confirm Import');
%    switch choice
%        case {'No','Cancel'}
%            fprintf(1,'TTTR Import canceled by user.\n');
%            return;
%    end
%end

%% Read files and import data

% TTTRdata(nfiles) = struct('Name',[],'HeaderPars',[],'Header',[],'Events',[]);

for ff=1:nfiles
    fid = fopen(filenames{ff});
    
    Magic = fread(fid, 8, '*char')';
    Magic = Magic(Magic~=0); % remove empty characters.
    if not(strcmp(Magic,'PQTTTR'))
        error('Magic code invalid, this is not a PTU file.');
    end
    Version = fread(fid, 8, '*char')';
    Version = Version(Version ~=0); % remove empty characters
%     fprintf(1,'Tag Version: %s\n', Version);
    
    % Store these values
    [~,TTTRdata(ff).Name] = fileparts(filenames{ff});
    TTTRdata(ff).MagicNumber = Magic;
    TTTRdata(ff).Version = Version;
    
    ReadingHeader = true;
    while ReadingHeader
        % read Tag Head
        TagName = fread(fid, 32, '*char'); %Name of Tag parameter
        TagName = (TagName(TagName ~= 0))'; % remove #0 and more more readable      
        TagIdx = fread(fid, 1, 'int32');    % Index for array elements (-1 for scalars)
        TagType = fread(fid, 1, 'uint32');   % Type flag determines how to read value below
                                            
%         if TagIdx > -1
%           EvalName = [TagName '(' int2str(TagIdx + 1) ')'];
%         else
%           EvalName = TagName;
%         end
%         fprintf(1,'\n   %-40s', EvalName);

        if strcmp(TagName, 'Header_End')
            ReadingHeader = false; % still read value below, but then will exit loop
        end

        % Retrieve value based on type of parameter
        switch TagType
            case tyEmpty8
                fread(fid, 1, 'int64'); % Carries no meaningful data, just to separate groupings. Ignore
                continue; % do not store any values for this type -- pass to next iteration
%                 fprintf(1,'<Empty>');
            case tyBool8
                TagInt = fread(fid, 1, '*int64');
                Value = logical(TagInt); 
%                 if Value
%                     fprintf(1,'TRUE');
%                 else
%                     fprintf(1,'FALSE');
%                 end
            case tyInt8
                Value = fread(fid, 1, '*int64'); % integer value
%                 fprintf(1,'%d', Value);
            case tyBitSet64
                Value = fread(fid, 1, '*uint64'); % unsigned integer value
%                 fprintf(1,'%X', Value);
            case tyColor8
                Value = fread(fid, 1, '*uint64'); % unsigned integer value
%                 fprintf(1,'%X', Value);
            case tyFloat8
                Value = fread(fid, 1, 'double');
%                 fprintf(1, '%e', Value);
            case tyFloat8Array
                TagInt = fread(fid, 1, '*int64'); % gives length of array (in bytes)
                Value = fread(fid, [TagInt/8,1], 'double'); % read array
%                 fprintf(1,'<Float array with %d Entries>', TagInt / 8);
%                 fseek(fid, TagInt, 'cof');
            case tyTDateTime
                TagFloat = fread(fid, 1, 'double');
                Value = datestr(datenum(1899,12,30)+TagFloat); % convert to date string
%                 Value = datenum(1899,12,30)+TagFloat; % keep as datenum
%                 fprintf(1, '%s', datestr(datenum(1899,12,30)+TagFloat)); % display as Matlab Date String
            case tyAnsiString
                TagInt = fread(fid, 1, '*int64');
                TagString = fread(fid, TagInt, '*char');
                TagString = (TagString(TagString ~= 0))';
                Value = TagString;
%                 fprintf(1, '%s', TagString);
            case tyWideString
                % Matlab does not support Widestrings at all, just read bytes and
                % remove the 0's (up to current (2012))
                TagInt = fread(fid, 1, '*int64');
                TagString = fread(fid, TagInt, '*char');
                TagString = (TagString(TagString ~= 0))';
                Value = TagString;
%                 fprintf(1, '%s', TagString);
            case tyBinaryBlob
                TagInt = fread(fid, 1, '*int64');
                Value = fread(fid,TagInt,'*uint8'); % read in TagInt bytes
%                 fprintf(1,'<Binary Blob with %d Bytes>', TagInt);
%                 fseek(fid, TagInt, 'cof');
            otherwise
                error('Illegal Type identifier found! Broken file?');
        end
        
        % Set field name and value in Tags structure
        if TagIdx>-1 % deal with array elements
            if ischar(Value)
                Tags.(TagName){TagIdx+1} = Value;
            else
                Tags.(TagName)(TagIdx+1) = Value;
            end
        else
            Tags.(TagName) = Value;
        end
        
    end
    
    TTTRdata(ff).Tags = Tags;    
    fprintf(1,'Header loaded successfully for file %u of %u.\n',ff,nfiles);
    % Now read in data
    fprintf(1,'Loading data...');
    Nevents = Tags.TTResult_NumberOfRecords;
    TTTRdata(ff).Events = fread(fid,Nevents,'*ubit32');
    fprintf(1,' complete!\n');
    
    fclose(fid);

end