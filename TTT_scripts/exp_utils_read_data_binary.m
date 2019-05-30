function r = exp_utils_read_data_binary(filename,varargin)
%
% exp_utils_read_binary(filename)
%
% Read data from binary file and organise into a struct of arrays
%

fid = fopen(filename,'r'); % Open file to read only

%%%%
%%%% Header info
%%%%
n = fread(fid,1,'uint8');          % # characters in keyname
r.keyname = fread(fid,n,'*char')';  % Read experiment keyname
n = fread(fid,1,'uint8');          % # characters in experiment name
r.exp_name = fread(fid,n,'*char')'; % Read experiment name
n = fread(fid,1,'uint8');          % # characters in experiment date
r.date = fread(fid,n,'*char')';     % Read experiment date
n = fread(fid,1,'uint8');          % # characters in experiment time
r.time = fread(fid,n,'*char')';     % Read experiment time
r.nrecords = fread(fid,1,'uint32');  % total # of records
r.nblocks = fread(fid,1,'uint32');   % total # of blocks
r.ntrials = fread(fid,1,'uint32');   % total # of trials

alldata = 1; % Flag to read all data
rblocks = 1:r.nblocks; % Vector of blocks to read
rtrials = 1:r.ntrials; % Vector of trials to read

%%%%
%%%% Parse optional arguments
%%%%
if nargin>1
  j = 1;
  while j<=(nargin-1)
    if strcmp(varargin{j},'blocks') % vector of blocks to read
      alldata = 0;
      rblocks = varargin{j+1};
      j = j + 2;
    elseif strcmp(varargin{j},'trials') % vector of trials to read
      alldata = 0;
      rtrials = varargin{j+1};
      j = j + 2;
    elseif strcmp(varargin{j},'alldata') % flag to read all data or only automated data
      alldata = varargin{j+1};           % 0 - only auto data; 1 - all data, manual + auto (default)
      j = j + 2;
    end;
  end;
end;

%%%%
%%%% Read in data
%%%%
dnlist = cell(100,1); % Keep record of datanames so that existing fields are not overwritten
dncount = 1; % Next available cell in dnlist
for j=1:r.nrecords  
  block = fread(fid,1,'uint32'); % Block #
  trial = fread(fid,1,'uint32'); % Trial #
  nbytes = fread(fid,1,'uint32'); % # bytes in record
  % Determine whether to read record
  if alldata
    read = 1;
  elseif (sum(rblocks==block) && sum(rtrials==trial))
    read = 1;
  else
    read = 0;
  end;
  % Read or skip record
  if read
    nentries = fread(fid,1,'uint32'); % # entries in record
    for k=1:nentries
      n = fread(fid,1,'uint8'); % # characters in dataname
      dataname = fread(fid,n,'*char')'; % Read dataname
      existflag = 0; % Flags whether dataname is in dnlist
      % Check whether dataname exists
      q = 1;
      while ((q<dncount) && (~existflag))
        if strcmp(dataname,dnlist{q})
          existflag = 1; % Dataname does already exist
        end;
        q = q + 1;
      end;
      % If dataname does not yet exist, create its cell array, and enter
      % dataname into dnlist
      if ~existflag
        dnlist{dncount} = dataname; dncount = dncount + 1;
        if alldata
          eval(['r.' dataname ' = cell(r.nrecords,1);']);
        else
          eval(['r.' dataname ' = cell(length(rblocks),length(rtrials));'])
        end;
      end;
      % Extract data and place it in the cell array for dataname
      n = fread(fid,1,'uint8'); % # characters in classname
      classname = fread(fid,n,'*char')'; % Read classname
      if strcmp(classname,'char'), classname = '*char'; end; % To read chars in as strings
      n = fread(fid,1,'uint8'); % # dimensions in data
      datasize = fread(fid,n,'uint32')'; % Size of data
      data = fread(fid,prod(datasize),classname); % Read in data
      data = reshape(data,datasize); % Reshape data into array of size datasize
      if alldata
        command = ['r.' dataname '{' num2str(j) '} = data;'];
      else
        ind1 = num2str(find(rblocks==block)); % Index into data array for this block
        ind2 = num2str(find(rtrials==trial)); % Index into data array for this trial
        command = ['r.' dataname '{' ind1 ',' ind2 '} = data;'];
      end;
      eval(command); % Enter data into struct
    end;
  else
    fseek(fid,nbytes,'cof'); % Move forward nbytes to next record
  end;  
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fwrite(fid,block,'uint32'); % Block #
% fwrite(fid,trial,'uint32'); % Trial #
% fwrite(fid,n,'uint32'); % # of remaining bytes in record
% fwrite(fid,nentries,'uint32'); % # of entries in record
% for j=1:nentries
%   fwrite(fid,length(dataname{j}),'uint8'); % # characters in dataname
%   fwrite(fid,dataname{j},'char'); % Write dataname  
%   fwrite(fid,length(classname{j}),'uint8'); % # characters in classname
%   fwrite(fid,classname{j},'char'); % Write classname
%   fwrite(fid,length(datasize{j}),'uint8'); % # of dimensions to data
%   fwrite(fid,datasize{j},'uint32'); % Size of each dimension
%   fwrite(fid,dataentry{j},classname{j}); % Enter data
% end;