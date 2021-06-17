% name: readsb
% author: Chris Proctor, SSAI / NASA GSFC OBPG
% syntax: [data, sbHeader, headerArray] = readsb(fileName, <optional argument pairings>)
% 
% readsb is designed to open and read data files that are in a SeaBASS 
% format (https://seabass.gsfc.nasa.gov/). Data outputs can be returned as either
% a cell array or as a structure using the names of the data field headers 
% (e.g. dataStruc.depth, dataStruc.chl, dataStruc.lw412, etc). Additional 
% (optional) outputs include a structure of the header, and cell arrays of 
% the header and/or data matrix.
%
% Inputs Arguments:
%   
%   FILENAME
%       Required.
%       Syntax: readsb(FILENAME)
%       FILENAME must be the first argument. It is the name of the single
%       SeaBASS-formatted file you wish to load.
%   
%   SetMissingValue  <false OR a non-zero number>
%       Optional.
%       Example: readsb(FILENAME, 'SetMissingValue', false)
%             or readsb(FILENAME, 'SetMissingValue', -9999);
%       By default, readsb converts missing values (as specified in the metadata
%       header, e.g., /missing=-9999) to NaN. Set this to false to preserve the
%       original number used for missing values, or set to a number of your
%       choice to replace all missing values in the data with your new number.
%   
%   SetBDLValue  <NaN or number>
%       Optional.
%       Example: readsb(FILENAME, 'SetBDLValue', NaN)
%             or readsb(FILENAME, 'SetBDLValue', -8888);
%       Set this to a numeric value (including NaN) to replace any "below_detection_limit" 
%       values ("BDL" values for short) in the data. 
%       BDL values are used as replacement values when a parameter was measured 
%       for but not found (i.e., concentrations were lower than could be detected 
%       by the instrument used.) BDL values are distinct from "missing" 
%       values as missing indicates that something wasn't measured or reported.
%
%   SetADLValue  <NaN or number>
%       Optional.
%       Example: readsb(FILENAME, 'SetADLValue', NaN)
%             or readsb(FILENAME, 'SetADLValue', -7777);
%       Set this to a numeric value (including NaN) to replace any "above_detection_limit" 
%       values ("ADL" values for short) in the data. 
%       ADL values are used as replacement values when a parameter was measured 
%       for but not found (i.e., concentrations were higher than the range
%       allowed by the instrument used.) ADL values are distinct from "missing" 
%       values as missing indicates that something wasn't measured or reported.
%
%   MakeStructure <TRUE>
%       Optional, FALSE by default.
%       Example: readsb(FILENAME, 'MakeStructure', true)
%       By default readsb outputs all values in a cell array. Set 
%       MakeStructure to true if you instead want data output as a structure. 
%       The structure will contain a field named for each data field 
%       (e.g. data.chl, data.phaeo, data.LU410, etc.) 
%       
%       Also, a special field called "wavelength_lists" is made
%       containing a structure of all wavelength-associated variable names each 
%       with a numerical array of wavelengths associated
%       with them. (e.g., 'Rrs' [410 443 480 555 670])
%
%   NewTextFields <single dimensional cell array of strings>
%       Optional.
%       Example: myNewFields = {'shipName', 'example2', 'example3'};
%                readsb(FILENAME, 'NewTextFields', myNewFields);
%             or readsb(FILENAME, 'NewTextFields', {'all'});
%       This argument is only used to optimize the reading speed of files
%       containing uncommon-to-SeaBASS columns of text data, if you can 
%       pre-identify those field names (otherwise the script will detect them for you). 
%       It is unnecessary to list commonly used text fields like "station",
%       "cruise", or "type", as those are already in the list of
%       defaultTextFieldNames located further below.
%
%       Entering just {'all'} as the input for NewTextFields will
%       cause the entire data block to be read as text. An example of why
%       you might use this option: if you just wanted to rearrange
%       or make minor edits to the file (doing little or no math), having all 
%       the columns as text simplifies the code needed to re-write it.
%   
%   ErrorIfDataUseWarning <TRUE>
%       Optional, FALSE by default.
%       Example: readsb(FILENAME, 'ErrorIfDataUseWarning', true);
%       If set to true, an error message will occur if the SeaBASS file contains the
%       data_use_warning header which contains a list of one or more special
%       warnings. Consult the metadata comments or accompanying docs for more info.
%       Such files may fine to use, but it depends on your application.
%       Standard data_use_warnings include:
%       - optically_shallow: the data may contain values impacted by
%                            bottom reflectance)
%       - experimental: the data contain non-geospatial measurements, such
%                       as incubations, laboratory measurements, or otherwise 
%                       manufactured/manipulated environmental conditions
%
%   QuietMode <TRUE>
%       Optional, FALSE by default
%       Example: readsb(FILENAME, 'QuietMode', true);
%       If set to true, certain output messages and warnings won't be
%       printed to the console. Thus, use with caution.
%
%   DatenumTime <FALSE>
%       Optional, TRUE by default.
%       Example: readsb(FILENAME, 'DatenumTime', false);
%       If true, "time" (in the data block, if present) is converted from HH:MM:SS to MATLAB's 
%       datenum format. Additionally, a new column (or structure field) will be
%       added to the "data" output containing MATLAB's datenum time for every
%       measurement row in the original file (if necessary, that
%       information will be obtained from the header). Set to false to
%       disable these features.
%
%   FillAncillaryData <TRUE>
%       Optional, FALSE by default.
%       Example: readsb(FILENAME, 'FillAncillaryData', true);
%       Set to true to force readsb to create fields for lat, lon, depth, station, and datenum in the data
%       output, based on the metadata header values (i.e., SeaBASS files 
%       leave such info in the metadata header if it is the same for the
%       entire file, but this duplicates the header info on every row, to
%       simplify your other scripts).
%
% Outputs:
%
%   A variable number of outputs (1-4) are supported. Sequentially they are:
%
%   [data, sbHeader, headerArray, dataArray] = readsb(FILENAME)
%
%   data (1st output) is by default a cell matrix, but will be a structure if the
%       'MakeStructure' argument is used.
%
%   sbHeader (2nd output) is a structure of the SeaBASS file's header where the fields
%       are based on metadata header names (e.g. sbHeader.missing).
%
%   headerArray (3rd output) is a cell array containing the verbatim contents of the SeaBASS
%       header. It can be useful if you are plan to re-write a file.
%
%   dataArray (the 4th output) is a cell array containing the data block.
%       Useful if you plan to minimally change the data and re-write
%       a SeaBASS file.
%
% Examples:
%   
%   [data, header] = readsb('myfile.csv');
%   [data, header, headerArray, dataArray] = readsb('myfile.txt', 'MakeStructure', true, 'FillAncillaryData', true, 'SetBDLValue', NaN, 'SetADLValue', NaN, 'ErrorIfDataUseWarning', true);
%
% Notes:
%
% * This function is designed to work with files that have been properly
% formatted according to SeaBASS guidelines (i.e. Files that passed FCHECK).
% Some error checking is performed, but improperly formatted input files
% could cause this script to error or behave unexpectedly. Files
% downloaded from the SeaBASS database should already be properly formatted, 
% however, please email seabass@seabass.gsfc.nasa.gov and/or the contact listed
% in the metadata header if you identify problems with specific files.
%
% * It is always HIGHLY recommended that you check for and read any metadata
%   header comments and/or documentation accompanying data files. Information 
%   from those sources could impact your analysis.
%
% * MATLAB compatibility: the latest readsb was tested with R2016b
%
% /*=====================================================================*/
%                  NASA Goddard Space Flight Center (GSFC) 
%          Software distribution policy for Public Domain Software
% 
% The readsb code is in the public domain, available without fee for 
% educational, research, non-commercial and commercial purposes. Users may 
% distribute this code to third parties provided that this statement appears
% on all copies and that no charge is made for such copies.
%
% NASA GSFC MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THE SOFTWARE
% FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
% WARRANTY. NEITHER NASA GSFC NOR THE U.S. GOVERNMENT SHALL BE LIABLE FOR
% ANY DAMAGE SUFFERED BY THE USER OF THIS SOFTWARE.
% /*=====================================================================*/

% Here is a simplified example of how to read through a large number of files
% arranged in a bundle downloaded using the SeaBASS File Search. When performing your 
% analysis, remember to also look at the relevant documents and metadata header comments.
%
% Start in the base folder and adapt further as needed: 
%
% fileList = dir('**/*');
% for ix = numel(fileList):-1:1
%    
%    % skip documents and directories since readsb.m isn't designed to open those
%    if contains(fileList(ix).folder, '/documents', 'Ignorecase', true) == true || (fileList(ix).isdir == true)
%        fileList(ix) = [];
%    % read data files
%    else
%        [data, header] = readsb([fileList(ix).folder, '/', fileList(ix).name], 'MakeStructure', true, 'FillAncillaryData', true, 'SetBDLValue', NaN, 'ErrorIfDataUseWarning', true);
%        % Do XYZ analysis 
%    end
%     
% end

% Changelog:
% 2013-03-19: Initial Release
% 2013-04-08: Fixed problem caused when reading certain files with columns of text data.
% 2013-06-10: Major revisions were made to this function to enable it to automatically read files containing non-standard fields of text data without requiring additional user input.
% 2014-06-04: Added support for files with 'below_detection_limit header and a new error for files that have no rows of data.
% 2015-11-17: Added "all" option for newTextFields to read all data as text.
% 2015-12-02: Added 'fillancillarydata' option. Set various useful options enabled by default. Refined input syntax to accept logicals.
% 2016-02-10: Fixed bug in 'fillancillarydata' option that could incorrectly set 'lon' & 'depth' values
% 2016-05-20: Added ADL functionality to mimic BDL functionality. Updated suggested nominal values for these to -7777 and -8888, respectively.
% 2017-04-11: Added support for reading "date_time" (contained only in SeaBASS "validation" files) and other minor updates
% 2017-08-31: Added 'readsb:TooManyOrFewOutputVars' files containing the header "/optical_depth_warning=true"
% 2017-10-26: MakeStructure now includes a new field called wavelength_lists(*changed, see later updates), which has a structure containing names of all products with wavelengths along with numerical arrays of those wavelengths
% 2018-01-10: Miscellaneous minor bug fixes and improvements (particularly related to detecting and replacing missing/BDL/ADL values). Swapped X & Y dimensions for the optional headerArray output. Also added an 
%             example of how to parse through a directory tree of SeaBASS files as is downloaded from the website's File Search, see it just above this changelog. 
% 2018-03-22: Added station to the values filled by the optional FillAncillaryData option
% 2018-04-02: Moved "wavelength_lists" from the data output to the header output, and it's now included for all "MakeStructure" choices. Added a header field 
%             called "fields_list" as a cell structure with each field name split into its own cell. Order is preserved relative to the column number in the original file, 
%             i.e., fields_list{1} = name of the contents of the 1st column
% 2019-04-15  Added support for times that contain fractional seconds
% 2019-07-24  Added support for /data_use_warning (to replace and expand upon /optical_depth_warning)
function varargout = readsb(fileName,varargin)

% Optional user argument/input validation and parsing
p = inputParser;
p.addRequired('fileName',@(x) ischar(x));

% These optional user input values must be true or false
validationFcn = @(x) isscalar(x) && (islogical(x) || x == 1 || x == 0);
p.addParameter('MakeStructure', false, validationFcn); % addParamValue changed by jscott on 2016/04/12 to addParameter to support R2016a
p.addParameter('DatenumTime', true, validationFcn);
p.addParameter('FillAncillaryData', false, validationFcn);
p.addParameter('ErrorIfDataUseWarning', false, validationFcn);
p.addParameter('QuietMode', false, validationFcn); % suppresses some messages from being printed to the console

% SetMissingValue user input value must be numeric. Disable a given option by providing
% the argument false or 0 (Zero). As an intended side effect, the replacement value
% can't be set to zero
validationFcn = @(x) (isnumeric(x) && isscalar(x)) || x == false;
p.addParameter('SetMissingValue', NaN, validationFcn);

% SetBDLValue user input must be numeric (including NaN)
validationFcn = @(x) isnumeric(x) && isscalar(x);
p.addParameter('SetBDLValue', false, validationFcn);

% SetADLValue user input must be numeric (including NaN)
validationFcn = @(x) isnumeric(x) && isscalar(x);
p.addParameter('SetADLValue', false, validationFcn);

% NewTextFields must be entered as a single dimensional cell array of characters
validationFcn = @(x) iscell(x) && all(cellfun(@ischar, x)) && min(size(x)) == 1;
p.addParameter('NewTextFields', false, validationFcn);

p.parse(fileName,varargin{:});
optionalArgs = p.Results;

% Following is a list of standardized field names whose data columns should 
% include text. This isn't an exhaustive list, though it is optional
% because readsb will automatically detect text columns. However, it can
% reduce run-time to pre-identify them either by permanently adding their
% names to this list, or using the 'NewTextFields' optional input argument.
knownTextFieldNames = {
     'associated_files','associated_file_types','cruise','chors_id','date_time',...
     'file','hpl_id','hplc_gsfc_id','instrument','sample','SN',...
     'station','time','type',...
     };
 
%% Read data

if optionalArgs.QuietMode == false
    disp(['Reading file: ', fileName]);
end

if ~exist(fileName, 'file')
   
    error('readsb:FileNotFound', 'Could not find or open file %s', fileName);

elseif exist(fileName, 'dir')
    
    error('readsb:FolderNotFile', 'Could not open file (%s). It is a directory', fileName);
      
end

% Read in the entire data block as text
fileID = fopen(fileName, 'rt');
    allDataLines = textscan(fileID, '%s', 'Delimiter', '\n'); % , 'BufSize', 131072); % Commented out by jscott on 2016/04/12 for compatibility with R2016a
fclose(fileID);

% Read meta-data headers
[header, headerArray, nHeaderLines] = readsbheader(allDataLines);
allDataLines = allDataLines{1}(nHeaderLines + 1:end);
nDataRows = length(allDataLines);

if nDataRows == 0
   error('readsb:NoData', 'No data rows found in the file %s', fileName); 
end

fields = regexp(header.fields, ',', 'split');
% Error if any fields (names) are empty. Ideally this never occurs in a valid file.
if any( cellfun('isempty', fields) )
    error('readsb:EmptyDataFieldName', ['Empty name found when reading "fields" ',...
          'metadata header: %s'], fileName);
end

% Obtain the missing value from the header. There should be only 1 missing
% value in normal files, but error checking is performed just in case there 
% are unusual files with multiple values.
missing = header.missing;

if ischar(missing)

    splitMissing = regexp(missing, ',', 'split');

    for ix = 1:length(splitMissing)
    
        if ~isempty( str2num(splitMissing{ix}) ) %#ok<ST2NM>
            multipleMissing{ix} = str2num(splitMissing{ix}); %#ok<AGROW,ST2NM>
        end
        
    end
    
    missing = splitMissing{1};
    
elseif size(missing,2) > 1

    % Create a "multipleMissing" cell array
    for ix = 1:length(missing)
        multipleMissing{ix} = missing(ix); %#ok<AGROW>
    end

    missing = missing(1);
     
end

if exist('multipleMissing', 'var')
    
    warning('MATLAB:readsb:MultipleMissingValues',['Multiple missing values found ',...
            'in the meta-data header. SeaBASS headers normally have only a ',...
            'single missing value.']);
    
end

% textMissing is a text-version of the missing values, needed later in the
% replaceMissing function to avoid rare cases when num2str causes a loss in
% numerical precision (e.g., rounding -99.99999 to -100)
textMissing = regexprep( headerArray{find( cellfun(@(x) ~isempty(x), regexp(headerArray,'^/missing=(.*)', 'match')) , 1 )}, '/missing=', '' );
textMissing = regexp(textMissing, ',', 'split');

% Interpret delimiter from header delimiter field (space, comma or tab)
whitespace = [' ' '\b' '\t'];
switch lower(strtrim(header.delimiter))
    case 'comma'
        delimiter = ',';
    case 'space'
        delimiter = whitespace;
    case 'tab'
        delimiter = whitespace;
    otherwise
        error('readsb:InvalidDelimiter','File delimiter not recognized (acceptable delimiters are: space, comma, or tab');
end

% add any user-defined text fields to the standard list
if iscell(optionalArgs.NewTextFields)
    
    % check for 'all' newTextFields option in which case the all data
    % fields will be interpreted as text
    if (numel(optionalArgs.NewTextFields) == 1) && (strcmpi(optionalArgs.NewTextFields,'all'))
        knownTextFieldNames = fields;
    else
        knownTextFieldNames(end+1:end+numel(optionalArgs.NewTextFields)) = optionalArgs.NewTextFields;% = [knownTextFieldNames, optionalArgs.NewTextFields];
    end
    
end

% find the column index of any known text fields in preparation for textscan
textColIndexes = getindex(fields, knownTextFieldNames);


% This loop exists so that if issues occur during textscan, each line can 
% be parsed to find the problematic (i.e. non-numeric) fields, allowing an
% informative and specific message about the issue to be printed for troubleshooting
warning('ON', 'MATLAB:readsb:TextDataDetected');
nLine = 1;
while nLine > 0 && nLine <= nDataRows

    [dataMatrix, textColIndexes, nLine] = readdata(fileName, nLine, fields, textColIndexes, nHeaderLines, delimiter, missing, allDataLines, optionalArgs.QuietMode);

end

%%
% Replacement (i.e., "missing") values are automatically replaced with NaN 
% (unless the user disables or specifies some other replacement number). 
% The SetBDLValue may be used to replace below detection limit values with 
% NaNs or other user specified values.
setMissing = true;
if optionalArgs.SetMissingValue == 0
    setMissing = false;
% warn user if zero happens to be the native missing value    
elseif missing == 0
    warning('MATLAB:readsb:badreplacementvalue','This file''s native "/missing" value is zero, but using zero is deprecated, and could cause unexpected issues, such as replacing valid zero data (e.g., surface measurement depth value at the surface))');
end


fileUsesBDL = find( cellfun(@(x) ~isempty(x), regexp(headerArray,'/below_detection_limit=.*', 'match')), 1 ); % check if file contains BDL, otherwise no need to attempt replacing
if ~isempty(fileUsesBDL) && ~islogical(optionalArgs.SetBDLValue) && optionalArgs.SetBDLValue ~= header.below_detection_limit
   
   setBDL = true;
   if header.below_detection_limit == 0
        warning('MATLAB:readsb:badreplacementvalue','This file''s native "/below_detection_limit" value is zero, but using zero is deprecated, and could cause unexpected issues, such as replacing valid zero data (e.g., measurement depth value at the surface))');
   end
   
else
    setBDL = false;
end

fileUsesADL = find( cellfun(@(x) ~isempty(x), regexp(headerArray,'/above_detection_limit=.*', 'match')), 1 ); % check if file contains ADL, otherwise no need to attempt replacing
if ~isempty(fileUsesADL) && ~islogical(optionalArgs.SetADLValue) && optionalArgs.SetADLValue ~= header.above_detection_limit
   
   setADL = true;
   if header.below_detection_limit == 0
        warning('MATLAB:readsb:badreplacementvalue','This file''s native "/above_detection_limit" value is zero, but using zero is deprecated, and could cause unexpected issues, such as replacing valid zero data (e.g., measurement depth value at the surface))');
   end
   
else
    setADL = false;
end

if (setMissing || setBDL || setADL)
    
    for nField = 1:length(dataMatrix)

        if iscell(dataMatrix{nField})
           isText = true;
        else 
           isText = false;
        end       
        
        if setMissing
            
            % Replace all missing values (there should be only 1 missing value
            % but this loop deals with multiple missing... just in case.)
            numMissingValues = 1;
            if exist('multipleMissing', 'var')
                numMissingValues = length(multipleMissing);
            else % in almost all cases there is just one missing value
                multipleMissing{1} = missing;
            end

            for numMissing = 1:numMissingValues

                if isText % if a column is text data, perform a direct string comparison to avoid rare cases of precision loss involving num2str conversion
                   multipleMissing{numMissing} = textMissing{numMissing};
                elseif ischar(multipleMissing{numMissing}) % for numeric data columns, ensure that the missing value has been properly parsed as numeric (required)
                    multipleMissing{numMissing} = str2double(multipleMissing{numMissing});
                end
                
                dataMatrix{nField} = replacemissingvalues(dataMatrix{nField}, isText, multipleMissing{numMissing}, optionalArgs.SetMissingValue);
                
            end
            
        end
        
        if setBDL
            
            if isText % get a text version of the value for search and replace
            	bdlValue = regexprep( headerArray{find( cellfun(@(x) ~isempty(x), regexp(headerArray,'/below_detection_limit=(.*)', 'match')) , 1 )}, '/below_detection_limit=', '' );                    
            else % for numeric data columns, ensure that the value has been properly parsed as numeric (required)
            	bdlValue = header.below_detection_limit;
            end
            
            dataMatrix{nField} = replacemissingvalues(dataMatrix{nField}, isText, bdlValue, optionalArgs.SetBDLValue);
        
        end
        
        if setADL
            
            if isText % get a text version of the value for search and replace
            	adlValue = regexprep( headerArray{find( cellfun(@(x) ~isempty(x), regexp(headerArray,'/above_detection_limit=(.*)', 'match')) , 1 )}, '/above_detection_limit=', '' );                    
            else % for numeric data columns
                adlValue = header.above_detection_limit;
            end
            
            dataMatrix{nField} = replacemissingvalues(dataMatrix{nField}, isText, adlValue, optionalArgs.SetADLValue);
        
        end
        
    end
    
    missing = optionalArgs.SetMissingValue;
    
end

 
% If 'DatenumTime' set: if the time field exists, reformat "time" (HH:MM:SS)
%  with MATLAB datenums. Additionally, create new field called datenum.
if logical(optionalArgs.DatenumTime) == true
    
    % find time index and update that field
    if  ~isempty(getindex(fields, 'time'))
        dataMatrix(getindex(fields, 'time')) = { sbtime2datenum( dataMatrix{getindex(fields, 'time')}, missing ) }; 
    end
    
    % create new datenum data
    dataMatrix(end+1) = {makedatenum( dataMatrix, header, lower(fields) )};
    header.fields = [header.fields,',datenum'];
    fields{end+1} = 'datenum';

end


% FillAncillaryData
if logical(optionalArgs.FillAncillaryData) == true
    
    [newData, newFields] = fillancillarydata(dataMatrix,header);
    if ~isempty(newFields)
        
        dataMatrix = [dataMatrix,newData];
        header.fields = [header.fields,sprintf(',%s',newFields{:})];
        fields(end+1:end+numel(newFields)) = newFields;
        
    end
    
    % make datenum if not already created
    if ~ismember(fields,'datenum')
        
        dataMatrix(end+1) = {makedatenum( dataMatrix, header, lower(fields) )};
        header.fields = [header.fields,',datenum'];
        fields{end+1} = 'datenum';

    end
    
end


% Calculate extra information about field-names and wavelengths: make a new structure field 
% that contains a structure of all the fields that contain wavelengths (if any). The basename 
% is listed along with a numerical array of all the wavelengths. Fields with wavelengths are
% determined from pattern matching. This extra information is added to
% the header.
% NOTE: _ex### _em### fields aren't included, not sure yet how to list those
%     also, this logic fails to separate suffixes, e.g., ag400 vs ag400_sd
patternSeabassWavelengths = '^[\w*]+(?<!_ex|em)\d{3,4}(\.\d{1+})?';
ixFieldsWithWavelengths = find(~cellfun(@isempty, regexp(fields,patternSeabassWavelengths, 'match')) == 1);

% determine the list of field basename (i.e., "Rrs" given multiple Rrs<lambdas>
wlFieldsBaseNames = unique( regexprep(fields(ixFieldsWithWavelengths), '\d{3,4}(\.\d+)?.*', '') );

% Parse through the fields with wavelengths. Determine the numerical
% wavelengths and assign them into a structure
for ix = 1:numel(wlFieldsBaseNames)

    baseNameIndexes = find( strncmpi(fields(ixFieldsWithWavelengths), wlFieldsBaseNames{ix}, length(wlFieldsBaseNames{ix})) );

    % list the wavelengths associated with a given field
    % Note that this does not work with fields using _ex _em, e.g., FL-A_ex488_ex670
    fieldWavelengths = regexprep(fields(ixFieldsWithWavelengths(baseNameIndexes)), '^[\w*]+?(\d{3,4}(\.\d+)?).*', '$1');

    wavelength_lists.(makestructurename(wlFieldsBaseNames{ix})) = unique(str2double(fieldWavelengths));
    %wavelength_lists.(wlFieldsBaseNames{ix}) = str2double(fieldWavelengths);

end

if exist('wavelength_lists', 'var')
    header.wavelength_lists = wavelength_lists;
end

header.fields_list = regexp(header.fields,',','split'); 


% Output data as a structure, only if optArg was set
if logical(optionalArgs.MakeStructure) == false

    data = dataMatrix;
    
else
    
    data = makestructure(dataMatrix,fields);
    
end

% provide visible warning if file potentially contains optically shallow (i.e., bottom
% reflectance) data, as those can't be used normally for traditional
% validation and algorithm development
% [*** 2019-07-24 the optical_depth_warning header will be deprecated in
% the future and replaced entirely by the more general "data_use_warning"
% (i.e., discontinue use of the optical_depth_warning, instead use the keyword
% "optically_shallow" as part of /data_use_warning]
if isfield(header,'optical_depth_warning') && strcmpi(header.optical_depth_warning,'TRUE')
    
    % ErrorIfOpticallyShallow
    if logical(optionalArgs.ErrorIfDataUseWarning) == true
        
        error('readsb:ErrorIfDataUseWarning','The optical_depth_warning flag is true for this file, indicating that it contains measurements in optically shallow conditions (i.e., bottom reflectance). Remove (or set to false) readsb''s ErrorIfDataUseWarning argument if you want to override this error and read the file.');
        
    else
        
        warning('The optical_depth_warning flag is true for this file, indicating that it contains measurements in optically shallow conditions (i.e., bottom reflectance). Use with caution or exclude it if performing validation or algorithm development.');

    end    
end

% check for the presence of any formally-defined data_use_warning terms, which
% alert users to one or more data points that may have been collected under
% special conditions. Users should consult docs or user metadata comments to learn more and to decide 
% if it is appropriate to use such data for their specific application.
if isfield(header,'data_use_warning') && ~any(isnan(header.data_use_warning(:))) %this checks specifically for a single nan value

    data_use_warnings = lower(regexp(header.data_use_warning,',','split'));

    if logical(optionalArgs.ErrorIfDataUseWarning) == true

        error('readsb:ErrorIfDataUseWarning',['File reading aborted due to detecting data_use_warning flag. File contains data that are:\n %s\n'...
            'To override, set ErrorIfDataUseWarning to false'],...
            strjoin(data_use_warnings,'\n '));
    else

        warning('data_use_warning flag detected. File contains data that are:\n %s\n', strjoin(data_use_warnings,'\n '));

    end
end



% readsb allows 1-4 outputs.
varargout = cell(1,nargout);

if nargout == 1
    varargout{1} = data;
elseif nargout == 2
    varargout{1} = data;
    varargout{2} = header;
elseif nargout == 3
    varargout{1} = data;
    varargout{2} = header;
    varargout{3} = headerArray;
elseif nargout == 4
    varargout{1} = data;
    varargout{2} = header;
    varargout{3} = headerArray;
    varargout{4} = allDataLines;
else
    error('readsb:TooManyOrFewOutputVars','Function was called with an improper amount of output variables. Please specify between 1 and 4 output variables.');
end

end


%% -----------------------------------------------
% readsbheader
% Gather data from a SeaBASS meta-data header into a structure (and array)
% SeaBASS headers vary in length but always consist of the format:
% /begin_header
% /headera=XYZ
% /headerb=XYZ
% /...
% ! optional comments (anywhere in header)
% /end_header
function [sbHeader, headerArray, nHeaderLine] = readsbheader(allData)

    % Initialize SeaBASS header structure
    sbHeader = struct('investigators', NaN,...
    'affiliations', NaN,...
    'contact', NaN,...
    'experiment', NaN,...
    'cruise', NaN,...
    'station', NaN,...
    'documents', NaN,...
    'calibration_files', NaN,...
    'data_type', NaN,...
    'data_status', NaN,...
    'start_date', NaN,...
    'end_date', NaN,...
    'start_time', NaN,...
    'end_time', NaN,...
    'north_latitude', NaN,...
    'south_latitude', NaN,...
    'east_longitude', NaN,...
    'west_longitude', NaN,...
    'cloud_percent', NaN,...
    'wind_speed', NaN,...
    'wave_height', NaN,...
    'water_depth', NaN,...
    'measurement_depth', NaN,...
    'secchi_depth', NaN,...
    'delimiter', NaN,...
    'missing', -9999,...
    'below_detection_limit', -8888,... 
    'above_detection_limit', -7777,...
    'fields', NaN,...
    'units',NaN,...
    'comments','');
    
    maxBadLines = 5;
    badHeaderLines = 0;
    
    nHeaderLine = 1;
    headerLine = allData{1}{1};
    headerArray = cell(size(allData{1},1), 1);
    
    if ~strncmpi(headerLine, '/begin_header', 13)
        
        badHeaderLines = badHeaderLines + 1;
        warning('Malformed header or non-SeaBASS file: File does not start with "/begin_header"');
            
    end    
    
    while nHeaderLine <= length(allData{1})

        if badHeaderLines > maxBadLines
            error('readsb:MalformedSeaBASSFormat:MoreThanMaxBadlines',...
                ['Couldn''t read file. Either it isn''t a SeaBASS file or the format is malformed because',...
                ' more than %d lines in the header couldn''t be parsed correctly.'], maxBadLines);
        end

        headerArray(nHeaderLine) = {headerLine};

        % Record the header lines in a structure:    
        % Record comments ( i.e. beginning with "!" )
        if strncmp(headerLine, '!', 1)
            
            % Only record if a comment if text follows the "!"
            if length(headerLine) > 1

                if isempty(sbHeader.comments)
                    sbHeader.comments = {char(headerLine)};
                else
                    sbHeader.comments = [sbHeader.comments; {char(headerLine)}];
                end

            end

        % Record header lines ( i.e. those beginning with "/" )
        elseif strncmp(headerLine, '/', 1)

            % Split the header line using '=' as the delimiter
            tokens = regexp(headerLine, '=', 'split', 'once');

            if strncmpi(tokens(1), '/end_header', 11)

                break;

            elseif(size(tokens,2) == 2) 

                name = makestructurename( lower(char(tokens(1))) );
                % Remove any trailing brackets (e.g. 12:00:01[GMT] --> 12:00:01)
                value = regexprep(strtrim( tokens{2} ), '(.+)\[.+\]$', '$1');

                % It is unnecessary to warn the user if a header
                % VALUE is long (but still short enough it isn't truncated).  
                % The warning is then turned back on in case the header 
                % NAME is too long, because that can cause truncation)
                identifier = 'MATLAB:namelengthmaxexceeded';
                warning('off', identifier);
                numValue = str2double(value);
                warning('on', identifier);

                if isnan( numValue )
                    if strcmpi(strtrim(value),'NaN') || strcmpi(strtrim(value),'NA') % seabass headers sometimes have "na" or "nan" as values. An exception to this logic is caused if /fields starts with "na" like "nadir"
                        sbHeader.(name) = NaN;
                    else
                        sbHeader.(name) = value;
                    end
                else % convert value to a number if it isn't text
                    sbHeader.(name) = numValue;
                end

            end

        % Warn user if the header line begins with some other character 
        else

            warning('Header line # %d is malformed. SeaBASS headers are required to begin with "/" or "!", but this line was: " %s "', nHeaderLine, headerLine);
            badHeaderLines = badHeaderLines + 1;

        end

        % Read next header line
        nHeaderLine = nHeaderLine + 1;
        headerLine = allData{1}{nHeaderLine};

    end    
        
    % strink allHeaderLines by deleting the excess empty cells

    headerArray(cellfun(@(x) isempty(x), headerArray)) = [];
    
end


%% -----------------------------------------------
% readdata
% This function uses textscan to read the data block of a SeaBASS file.
% It returns the data matrix, as well as textColIndexes which is updated
% if any (unexpected) text data were encountered.
function [dataBlock, textColIndexes, nLine] = readdata(fileName, nLine, fields, textColIndexes, nHeaderLines, delimiter, missing, allDataLines, quietOption)
fileID = fopen(fileName, 'rt');

% specify formats used to read data
formatSpec = repmat({'%f'}, 1, length(fields));
formatSpec = sprintf('%s', formatSpec{:});

for columnPosition = 1:length(textColIndexes)
    formatSpec( textColIndexes(columnPosition) * 2-1:textColIndexes(columnPosition) * 2 ) = '%s';
end

if ~isempty(strfind(delimiter, ' '))
    regexpDelimiter = '\s+';
else
    regexpDelimiter = delimiter;
end

try
    % If "/missing=" value is a number, run textscan without TreatAsEmpty
    if ~ischar(missing)
        
        dataBlock = textscan(fileID, formatSpec, 'HeaderLines', nHeaderLines,...
        'Delimiter', delimiter, 'MultipleDelimsAsOne', 1, 'ReturnOnError', 0);
        
    else

        dataBlock = textscan(fileID, formatSpec, 'HeaderLines', nHeaderLines,...
        'Delimiter', delimiter, 'MultipleDelimsAsOne', 1, 'ReturnOnError', 0,...
        'TreatAsEmpty', {missing}); % Prior to MATLABR2016b "treatasempty" was not a cell array

    end
    
    % This section double checks the length of the data block. If textscan
    % expected a column of numbers but finds certain text symbols (like -),
    % (e.g. "123-456") it can silently overestimate or misinterpret the
    % number of rows, causing the data to become out of alignment. 
    % To address those cases, it is necessary to parse through the file a row 
    % at a time using regular expressions to detect the problem. This process
    % could theoretically be slow if the file is large (e.g. hyperspectral measurements).
    if length(dataBlock{1}) > length(allDataLines)
        
        if quietOption == false
            warning('MATLAB:readsb:TextDataDetected','Text data was detected that required extra parsing. No user action required, unless you wish to use the NewTextFields argument to speed future calls.');
        end
        % Disable this warning for the rest of the run; it is enabled earlier in
        % the run just in case it was disabled during the last run.
        warning('OFF', 'MATLAB:readsb:TextDataDetected'); 
        
        problemColumnLocated = false;

        while (problemColumnLocated == false) && (nLine <= length(allDataLines))
            
            dataLine = regexp(allDataLines{nLine}, regexpDelimiter, 'split');
            [problemColumnLocated, textColIndexes] = parserow(dataLine, nLine, missing, textColIndexes, fields);
            nLine = nLine + 1;        
            
        end

    % if lengths are equal, and no other errors, the file was read successfully. Hooray!
    else
        nLine = -1;
    end

% Text data in the data matrix is typically the cause of the following errors
catch errorMsg 
     
    % check if the error was caused by text data
    if ~isempty( strfind(errorMsg.message, 'Trouble reading floating point number from file (row') ) || ~isempty( strfind(errorMsg.message, 'Trouble reading ''Numeric'' field from file') )
        %badRow = str2double( regexp(errorMsg.message, '(?<=row )\d+','match') ); %this used to work in earlier matlab versions
        badRow = str2double( regexp(errorMsg.message, '\d+', 'match','once') ); % first number in the message is the row
    else
        % report other types of errors
        error(errorMsg.message);
    end
    
    % First check if badRow is in bounds. If it is, proceed and check if
    % textscan's error message correctly identified the row number
    % where the problem is (i.e. saving us the trouble of parsing for it)
    if badRow <= length(allDataLines)
    
        dataLine = regexp(allDataLines{badRow}, regexpDelimiter, 'split');
        [problemColumnLocated, textColIndexes] = parserow(dataLine, badRow, missing, textColIndexes, fields);    
        
        % However, if "badRow" had no problems (i.e. no unexpected text data), then
        % textscan was unsuccessful at pinpointing the row with bad data so
        % it is necessary to manually find it by going through the file 
        % line-by-line starting from the beginning (progress tracked by nLine) .
        if problemColumnLocated == false
            
            while (problemColumnLocated == false) && (nLine <= length(allDataLines))
                
                dataLine = regexp(allDataLines{nLine}, regexpDelimiter, 'split');
                [problemColumnLocated, textColIndexes] = parserow(dataLine, nLine, missing, textColIndexes, fields);
                nLine = nLine + 1;
                
            end
            
        end
    
    % If the badRow index was greater than the length of allDataLines then
    % that indicates the file contains text data AND there is a delimiting problem 
    % (in textscan) due to text characters like "-" or ":". The file must be parsed
    % line-by-line. Columns with values like 00-123 really toss a wrench in the works.
    else
        
        if quietOption == false
            warning('MATLAB:readsb:TextDataDetected','Text data was detected that required extra parsing. No user action required, unless you wish to use the NewTextFields argument to speed future calls.');
        end
        % Disable this warning for the rest of the run; it is enabled earlier in
        % this function in case it was disabled during the last run.
        warning('OFF', 'MATLAB:readsb:TextDataDetected'); 
        
        problemColumnLocated = false;

        while (problemColumnLocated == false) && (nLine <= length(allDataLines))
            
            dataLine = regexp(allDataLines{nLine}, regexpDelimiter, 'split');
            [problemColumnLocated, textColIndexes] = parserow(dataLine, nLine, missing, textColIndexes, fields);
            nLine = nLine + 1;

        end
        
    end
    
    dataBlock = -1;
    
end

fclose(fileID);

end


%% -----------------------------------------------
% getindex
% Pass in a list of the fields and the names of target fields and their
% indexes (i.e. numbered column positions) within the fields list are returned
function indexes = getindex(allFields, targetFields)

    indexes = find( ismember(lower(allFields), lower(targetFields)) );

end


%% -----------------------------------------------
% parserow
% This function will search a data row for any text data. It requires a row 
% of data, the index of that row (nLine), indexes of text columns and fields 
function [problemColumnLocated, textColIndexes] = parserow(dataLine, nLine, missing, textColIndexes, fields)

    problemColumnLocated = false;

    for nColumn = 1:length(dataLine)

        % check each value to see if it is 1) non-numeric, 2) not the missing value, and 3) not already a column that is recognized as text
        if isnan(str2double(dataLine(nColumn))) && ~strcmp(dataLine(nColumn), num2str(missing)) && ~ismember(nColumn, textColIndexes)

            fprintf('Non-numeric value (%s) detected in column %d (%s), row %d, so that column''s data was treated as text. No user action required unless you wish to use the NewTextFields argument to speed future calls.\n',dataLine{nColumn}, nColumn, fields{nColumn}, nLine );
            % Add the detected column to the fields list.
            textColIndexes = [textColIndexes nColumn]; %#ok<AGROW>
            problemColumnLocated = true;

        end

    end

end


%%  -----------------------------------------------
% sbtime2datenum
% convert the "time" column into a MATLAB datenum
function formattedTime = sbtime2datenum(time, missing)
    
    formattedTime = NaN(numel(time), 1);
    valid = not( strcmp(time,num2str(missing)) );
    formattedTime(valid,1) = datenum(time(valid,1));
    formattedTime(valid,1) = formattedTime(valid,1) - floor(formattedTime(valid,1));
    formattedTime(~valid,1) = missing;
    
end


%%  -----------------------------------------------
% makedatenum
% create new data columns containing the datetime converted to MATLAB datenum format
function newDatenum = makedatenum( dataMatrix, header, fields )
    
    numRows = numel(dataMatrix{1});
    dates = NaN(numRows, 1);
    times = NaN(numRows, 1);
    newDatenum = NaN(numRows, 1);
    
    % determine how to construct the datenum from the different fields that
    % can be included in a SeaBASS file. In fields, it consists of some combo of:
    % date time year hour day minute (second) headerstartdate headerstartime
    % Additionally, "date_time", not a true SeaBASS field, is supported because of 
    % its presence in SeaBASS validation files
    
    % date
    if any(strcmpi(fields, 'date'))  %(from data)
        
        try
            if iscell(dataMatrix{getindex(fields, 'date')})
                dformat = repmat({'yyyymmdd'},numRows,1);
                dindex = dataMatrix{getindex(fields, 'date')};
                dates(:,1) = cellfun(@datenum, dindex, dformat );
            else %this is most common
                dates(:,1) = datenum(num2str(dataMatrix{getindex(fields, 'date')}), 'yyyymmdd'); 
            end
            
        catch errorMsg
        end
        
    elseif sum(ismember(fields,{'year','month','day'})) == 3 % i.e., year & month & day are present
        
        try           
            if iscell(dataMatrix{getindex(fields, 'year')})
                iyear = getindex(fields, 'year');
                imonth = getindex(fields, 'month');
                iday = getindex(fields, 'day');
                for k = 1:numRows
                    dates(k,1) = datenum(sprintf('%4s%02s%02s', dataMatrix{iyear}{k},dataMatrix{imonth}{k},dataMatrix{iday}{k}),'yyyymmdd');
                end
            else
                dates(:,1) = datenum( dataMatrix{getindex(fields, 'year')}, dataMatrix{getindex(fields, 'month')}, dataMatrix{getindex(fields, 'day')} );
            end
        catch errorMsg
        end
    
    elseif any(strcmpi(fields, 'date_time'))
        
        try
            dformat = repmat({'yyyy-mm-dd HH:MM:SS'},numRows,1);
            dindex = dataMatrix{getindex(fields, 'date_time')};
            dates(:,1) = cellfun(@datenum, dindex(:,1), dformat(:,1) );
        catch errorMsg
        end
        
    else % get date from header
        
    	dates(:,1) = datenum( num2str(header.start_date), 'yyyymmdd');  

    end
    
    
    % time
    if any(strcmpi(fields, 'time'))

    	times(:,1) = dataMatrix{getindex(fields, 'time')};

    elseif sum(ismember(fields,{'hour','minute','second'})) >= 2  % hour&minute&(second) (from data fields)
        
        nzeros = zeros(numRows, 1);
        
        if any(ismember(fields,'second'))
            if ~iscell(dataMatrix{getindex(fields, 'second')})
                seconds = dataMatrix{getindex(fields, 'second')};
            else
                seconds = sprintf('%02d\n',dataMatrix{getindex(fields, 'seconds')}{:});
            end
        else %sometimes seconds-precision wasn't measured
            if ~iscell(dataMatrix{getindex(fields, 'hour')}) %see if hour was cell format and match/assume same format for seconds
                seconds = nzeros;
            else
                seconds = repmat({'00'},numRows,1);
            end
        end
        
        try
            if iscell(dataMatrix{getindex(fields, 'hour')}) && iscell(dataMatrix{getindex(fields, 'minute')})
                ihour = getindex(fields, 'hour');
                iminute = getindex(fields, 'minute');
                for k = 1:numRows
                    times(k,1) = rem( datenum(sprintf('%02s%02s%02s', dataMatrix{ihour}{k},dataMatrix{iminute}{k},seconds{k}),'HHMMSS'), 1);
                end
            else
                times(:,1) = datenum(nzeros, nzeros, nzeros, dataMatrix{getindex(fields, 'hour')}, dataMatrix{getindex(fields, 'minute')}, seconds);
            end
        catch errorMsg
        end
        
    elseif any(strcmpi(fields, 'date_time'))
        % date_time already has both date and time stored in "date"
        times(:,1) = 0;
        
    else % get time from header
        
        if isnan(header.start_time)
            warning('MATLAB:readsb:DatenumTime:MissingTimeInfo','Measurement time not found in data block nor start_time header. The new datenum data could be compromised, so check those values for accuracy before using them.');
            times(:,1) = 0;
        else
            times(:,1) = datenum( [0,0,0, str2double(strsplit(header.start_time,':'))]);
        end
    end
    
    if ~exist('errorMsg','var')
        newDatenum = dates + times;
    else
        warning('Errors occurred when trying to create a new "datenum" column :\n%s',errorMsg.message);
    end  
    
end


%%  -----------------------------------------------
% fillancillarydata
% Adds additional fields to the data block by finding values in the headers
% and duplicating them into their own columns (if not already present)
function [newData, newFields] = fillancillarydata(dataMatrix, header)
    
    % pairs of ("/fields" names : equivalent metadata-header-name)
    ancillaryFields = [{'lat','north_latitude'};...
                       {'lon','east_longitude'};...
                       {'depth','measurement_depth'};...
                       {'station','station'}];

    fields = lower(regexp(header.fields, ',', 'split'));
    ancFieldIndexes = ~ismember((ancillaryFields(:,1)),fields);
    numFields = sum(ancFieldIndexes);
    newFields = ancillaryFields(ancFieldIndexes,1);
    newData = repmat( {NaN(numel(dataMatrix{1}),1)},1,numFields );

    % create any needed new fields
    newFieldIndexes = find(ancFieldIndexes == 1); % find the indexes of ancfields that need to be created
    for ix = 1:numel(newFields)
        
        % fill station values as cell, else use doubles for the other numeric
        % metadata measurements
        if strcmp(ancillaryFields{newFieldIndexes(ix),2},'station') == true
            
            newData{ix}      = cell(numel(dataMatrix{1}), 1);
            newData{ix}(:,1) = {header.(ancillaryFields{newFieldIndexes(ix),2})};
            
        else
        
            newData{ix}(:,1) = header.(ancillaryFields{newFieldIndexes(ix),2});
            
        end
        
    end

end


%% -----------------------------------------------
% makestructurename
% Turns mixed text into a string appropriate for a structure's field name.
% Any non-word or non-numeric characters are replaced with an underscore
% and leading and trailing underscores are removed. 
% (e.g. _Custom-Field.Name+ becomes Custom_Field_Name)
function fieldName = makestructurename(value)

    pattern = '[^\w]';
    fieldName = regexprep(lower(value), '*', '_star_'); % added for rare variables with asterisks, like a*ph
    fieldName = regexprep(fieldName, pattern, '_');
    fieldName = regexprep(fieldName, '_+$', '');
    fieldName = regexprep(fieldName, '^_+', '');
    % the prefix 'n' is added any field names beginning with a digit
    fieldName = regexprep(fieldName, '^(\d)', 'n$1');

end


%% -----------------------------------------------
% makestructure
% Converts a cell array into a structure using the field names from the
% metadata header as the names of the structure's fields
function newStructure = makestructure(cells, fields)
newStructure = struct;
fields = makestructurename(fields);

numCols = length(fields);
if length(unique(fields)) < numCols
    warning('MATLAB:readsbheaders:duplicateFieldNames',['This file contains ',...
        'duplicate field names.\nConsequently, its structure field names will ',...
        'be renamed with incrementing numbers (using a suffix of _n#) ',...
        'to avoid overwriting the duplicated field names']);
end

% Assign data into named structure fields. A few SeaBASS files might
% contain duplicate field names, so additional instances of the same name
% are given an incrementing numerical suffix (e.g. chl, chl_n1, chl_n2). If
% the program detects more than maxTries duplicates of a name, it will end
% in an error (that situation is unexpected to ever occur in a normal
% SeaBASS file.)
for col = 1:numCols
    
    if ~isfield(newStructure, fields{col})
        
        newStructure.(fields{col}) = cells{col};
        
    else
        
        suffix = 2;
        maxTries = 10;
        while isfield(newStructure, [fields{col}, '_n', num2str(suffix)])
            
            suffix = suffix + 1;
            if suffix > maxTries
                  
                 error('readSB:MakeStructure:tooManyNameDuplicates',...
                     ['A field name is duplicated more than %d times in this file. ',...
                     'Please verify the /fields line in the file is correct'], maxTries);
                 
            end
            
        end
        
        newStructure.([fields{col}, '_n', num2str(suffix)]) = cells{col};

    end
    
end

end


%% ---------------------------------------------
% replacemissingvalues
% SeaBASS files contain a header field called /missing with a single value 
% that represents missing values in the file. A number like -9999 is often 
% used. The following function is provided to replace whatever the missing 
% values are with a numeric value of the users choice (NaN included).
% This function can also be used to make replacements for /below_detection_limit
function newCells = replacemissingvalues(cells, isText, oldMissingVal, newMissingVal)

if isText % If the data column is made up of text
    
	% error checking: assessing equality (needed for replacing numbers) is difficult to implement
    % for certain decimal precision, e.g., -99.99999 gets rounded to -100.
    % This prevents users from picking that option since the solutions are
    % too time consuming to be worth the effort.
    if ~isnan(newMissingVal) && ( newMissingVal ~= str2double(num2str(newMissingVal)) )
        disp(newMissingVal);
        error('readsb:SetMissingValue:BadReplaceValueChoice', 'Your SetMissingValue (%s) is a bad choice because its decimal precision causes equality assessment issues. Please pick an integer instead', newMissingVal);
    end

    % Note that this replacement has to account for some exception cases because some files have mismatches 
    % between trailing significant zeroes, i.e.,-999 needs to be treated as
    % the equivalent of -999.000 (i.e., via string/text comparison) 
    
    % Check for an integer missing value, in which case, allow excess trailing zeroes
    % in the data matrix to count as pattern matches
    if fix(str2double(oldMissingVal)) == str2double(oldMissingVal)

        pattern = regexprep(oldMissingVal,'^([-+]?\d+)(\.0*)?','^$1(\.0*)?$'); % allow trailing zeroes
        cells = regexprep(cells,pattern,num2str(newMissingVal));
        
    else

        cells = regexprep(cells,['^',oldMissingVal,'$'],num2str(newMissingVal), 'ignorecase');

    end
        
else % if the data column is numeric
    
    cells(cells == oldMissingVal) = newMissingVal;
    
    % textscan's treatasempty might have created NaN values, so if the
    % newMissingVal needs to be something other than NaN, so they are
    % automatically checked for in addition to the other oldMissingVal
    %if ~isnan(newMissingVal)
    %    cells(isnan(cells)) = newMissingVal;
    %end
  
end

newCells = cells;

end
