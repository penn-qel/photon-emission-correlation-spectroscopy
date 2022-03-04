function histoData = hist_import_PHU(filename)

% Read_PHU.m    Read PicoQuant Unified Histogram Files
% This is demo code. Use at your own risk. No warranties.
% Marcus Sackrow, Michael Wahl, PicoQUant GmbH, December 2013
% Updated by Michael Wahl, PicoQUant GmbH, June 2016

% Note that marker events have a lower time resolution and may therefore appear 
% in the file slightly out of order with respect to regular (photon) event records.
% This is by design. Markers are designed only for relatively coarse 
% synchronization requirements such as image scanning. 

% T Mode data are written to an output file [filename].out 
% We do not keep it in memory because of the huge amout of memory
% this would take in case of large files. Of course you can change this, 
% e.g. if your files are not too big. 
% Otherwise it is best process the data on the fly and keep only the results.

% All header data are introduced as Variable to Matlab and can directly be
% used for further analysis

% some constants
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
  
    %%
    
    fid=fopen(filename);
    
    fprintf(1,'\n');
    Magic = fread(fid, 8, '*char');
    if not(strcmp(Magic(Magic~=0)','PQHISTO'))
        error('Magic invalid, this is not a PHU file.');
    end
    Version = fread(fid, 8, '*char');
    Version = Version(Version ~= 0);
%     fprintf(1,'Tag Version: %s\n', Version);

    histoData.MagicNumber = Magic;
    histoData.Version = Version;
    % there is no repeat.. until (or do..while) construct in matlab so we use
    % while 1 ... if (expr) break; end; end;
    ReadingHeader = true;
%     counter = 0;
    while ReadingHeader
        % read Tag Head
        TagName = fread(fid, 32, '*char'); % TagHead.Ident
        TagName = (TagName(TagName ~= 0))';% remove #0 and make more readable
        TagIdx = fread(fid, 1, 'int32') ;   % TagHead.Idx
        TagType = fread(fid, 1, 'uint32');   % TagHead.Typ
                                            % TagHead.Value will be read in the
                                            % right type function  
%         if TagIdx > -1
%           EvalName = [TagIdent '(' int2str(TagIdx + 1) ')'];
%         else
%           EvalName = TagIdent;
%         end
%         fprintf(1,'\n   %-40s', EvalName);  
%         counter = counter +1
        if strcmp(TagName, 'Header_End')
%             disp('the end')
            ReadingHeader = false; % still read value below, but then will exit loop
        end

        % check Type of Header
        switch TagType
            case tyEmpty8
                fread(fid, 1, 'int64');   
%                 fprintf(1,'<Empty>');
            case tyBool8
                TagInt = fread(fid, 1, 'int64');
                Value = logical(TagInt); 
%                 if TagInt==0
%                     fprintf(1,'FALSE');
%                     eval([EvalName '=false;']);
%                 else
%                     fprintf(1,'TRUE');
%                     eval([EvalName '=true;']);
%                 end            
            case tyInt8
                Value = fread(fid, 1, 'int64');
%                 fprintf(1,'%d', Value);
            case tyBitSet64
                Value = fread(fid, 1, 'int64');
%                 fprintf(1,'%X', Value);
            case tyColor8    
                Value = fread(fid, 1, 'int64');
%                 fprintf(1,'%X', Value);
            case tyFloat8
                Value = fread(fid, 1, 'double');
%                 fprintf(1, '%e', Value);
            case tyFloat8Array
                TagInt = fread(fid, 1, 'int64'); % gives length of array (in bytes)
                Value = fread(fid, [TagInt/8,1], 'double'); % read array
%                 fprintf(1,'<Float array with %d Entries>', TagInt / 8);
%                 fseek(fid, TagInt, 'cof');
            case tyTDateTime
                TagFloat = fread(fid, 1, 'double');
                Value = datestr(datenum(1899,12,30)+TagFloat); % convert to date string
%                 Value = datenum(1899,12,30)+TagFloat; % keep as datenum
%                 fprintf(1, '%s', datestr(datenum(1899,12,30)+TagFloat)); % display as Matlab Date String
            case tyAnsiString
                TagInt = fread(fid, 1, 'int64');
                TagString = fread(fid, TagInt, '*char');
                TagString = (TagString(TagString ~= 0))';
                Value = TagString;
%                 fprintf(1, '%s', TagString);
            case tyWideString 
                % Matlab does not support Widestrings at all, just read and
                % remove the 0's (up to current (2012))
                TagInt = fread(fid, 1, 'int64');
                TagString = fread(fid, TagInt, '*char');
                TagString = (TagString(TagString ~= 0))';
                Value = TagString;
%                 fprintf(1, '%s', TagString);
            case tyBinaryBlob
                TagInt = fread(fid, 1, 'int64');
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

    
% Read all histograms into one matrix
for i = 1:Tags.HistoResult_NumberOfCurves
%     fseek(fid,Tags.HistResDscr_DataOffset(i),'bof');
    Counts(:,i) = fread(fid, Tags.HistResDscr_HistogramBins(i), 'uint32');
end    

histoData.rawCounts = Counts;
histoData.counts = nonzeros(Counts);
histoData.tags = Tags;

fclose(fid);
end