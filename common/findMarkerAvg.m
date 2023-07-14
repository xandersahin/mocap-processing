function [] = findMarkerAvg(trialNames,MARKERNAME)



%Select files to which the HJC locations should be added (multiple must be
%selected).  The same directory will be used to create a file for new files
%with HJC locations.
files = trialNames.trials;
directory = trialNames.basepath ;

% files=strvcat(files(:,:));


directory=strvcat(directory);
% infile=files(1,:); 
inpath=directory;


%Reformat file list and directory information

[a b]=size(trialNames.trials);

%read each *.trc file and append HJC markers
for j=1:b;
    infile=files{j};
    inpath=directory;
    x=0; mnames=0; tx=0; marks=0; time=0;
    [x,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(['\' infile '.trc'],inpath);   
    
    %compute time from the sampling frequency
    [a,b]=size(x);
    for i=1:a;
        time(i,1)=i/sfx-1/sfx;
    end
     
    inds.MARKER = (strmatch(MARKERNAME,mnames)-1)*3+1 ;

    MARKER = x(:,inds.MARKER:inds.MARKER+2);
    
    
    average = mean(MARKER);

    sprintf(['The average marker coordinates for ',files{j},' ',MARKERNAME,' is ', num2str(average), ' mm.'])
    

    
end

end %main function

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [pos,time,f,n,nmrk,mrk_names,file,inpath]=load_trc(infile,inpath)
%   [pos,time,f,n,nmrk,mrk_names]=load_trc(infile)
%   LOAD_TRC is used to open a data file from Motion Analysis Realtime
%   output (*.trc).
%
%   Inputs:
%       infile - trc file to be loaded
%                If infile is unspecified, the user is prompted to select 
%                the input file
%       inpath - directory of location where data file is located when no
%                path is specified, it defaults to current directory
%
%   Outputs:
%       pos - contains the meaured marker positions in order of the markers
%             that is columns 1-3 are the x,y,z components of marker 1
%                     columns 4-6 are the x,y,z components of marker 2
%                     ....
%       time - column vector of time
%       f - sample frequency
%       n - number of data frames
%       nmrk - number of markers
%       mrk_names - marker names
%
%   Updated: Feb. 15, 2006 (JWF)
%
%   MATLAB Version 7.1

n = nargin;
if (n==0);
    [infile, inpath]=uigetfile('*.trc','Select input file');
    if infile==0
        f='';
        n='';
        nmrk='';
        mrk_names='';
        data=[];
        return;
    end
    fid=fopen([inpath infile],'r');
    file = infile(1:length(infile)-4);
elseif (n==1);
    file = infile(1:length(infile)-4);
    fid=fopen(infile,'r');
else (n==1);
    file = infile(1:length(infile)-4);
    fid=fopen([inpath infile],'r');
end

if (fid==-1)
    disp('File not found');
    error(['We tried to find ' inpath infile ' but it wasnt in the directory. Check for it.'])
    f='';
    n='';
    nmrk='';
    mrk_names='';
    data=[];
    return;
end

% disp(['Lwroading file...' infile] );

%disregard header info
for h=1:2
    hdr=fgetl(fid);
end
file_info=fscanf(fid,'%f');
f=file_info(1);
nmrk=file_info(4);
hdr=fscanf(fid,'%s',4);
line=fgetl(fid);
line=fgetl(fid);
j=1;
jl=length(line);
for i=1:(nmrk+2)
    name=sscanf(line(j:jl),'%s',1);
    ii=findstr(line(j:jl),name);
    j=j+ii(1)+length(name);
    if i>2
        mrk_names(i-2,1)=cellstr(name);
    end
end

%mrk_names=fscanf(fid,'%s',nmrk+2);
for h=1:2
    hdr=fgetl(fid);
end

line=[];
data=[];

try
    while(length(data)<((nmrk*3)+2))
        line=fgetl(fid);
        data=sscanf(line,'%f');
    end
catch
    error(['There are probably extra unnamed markers in this file: ' inpath infile])
end

time(1,1)=data(2);
pos(1,:)=data(3:length(data));
i=1;
while feof(fid)==0
    i=i+1;
    line=fgetl(fid);
    data=sscanf(line,'%f');
    time(i,1)=data(2);
    for j=3:length(data)
        pos(i,j-2)=data(j);
        %%%%%%%%
    end
end

[n,nc]=size(pos);
if n==1
    time=1;
else
    time=time(1:n,1);
end
% Return the position data in m
% pos=pos/1000;

% Check if number of markers and data columns is the same

if size(pos,2) ~= nmrk
    pos = pos(:,1:nmrk*3) ;
%     disp('There were too many marker positions in file - likely virtual markers. Deleted extra entries.')
end

end

%-------------------------------------------------------------------------%

