function [] = offsetMarkers6mm(trialNames,outputPath)


%file and folder naming    
% userDir = getuserdir ;
% tempOutputPath = [userDir '\Temp\Sub' num2str(round(rand*1000000)) '\'] ;
% finalFolderName = '\TRCs_w_HJCs_offset';
% outputPathFull = [outputPath,finalFolderName]; 
outputPathFull = outputPath;

%Select files to which the HJC locations should be added (multiple must be
%selected).  The same directory will be used to create a file for new files
%with HJC locations.
files = trialNames.trials ;
directory = trialNames.basepath ;

% files=strvcat(files(:,:));

% if length(files)==1;
%     display('No Files Selected to offset markers.  Try Again.');
%     return;
% end

directory=strvcat(directory);
% infile=files(1,:); 
inpath=directory;


%Reformat file list and directory information

[a b]=size(trialNames.trials);
write_directory=outputPathFull;

%read each *.trc file and append HJC markers
for j=1:b;
    infile=files{j};
    inpath=directory;
    x=0; mnames=0; tx=0; marks=0; time=0;
    [x,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(['\' infile '_w_HJCs.trc'],inpath);   
    
    %compute time from the sampling frequency
    [a,b]=size(x);
    for i=1:a;
        time(i,1)=i/sfx-1/sfx;
    end
     
    % Add regression based HJCs

    mrknames= {'C7'; 'CLAV';
             'RSHO'; 'LSHO';
             'RUPA'; 'LUPA';
             'RFRM'; 'LFRM';
             'RELB'; 'LELB';
             'RWRA'; 'LWRA';
             'RWRB'; 'LWRB';
             'RHJC'; 'LHJC';
             'RASI'; 'LASI';
             'RPSI'; 'LPSI';
             'RTHI'; 'LTHI';
             'RKNE'; 'LKNE';
             'RMKN'; 'LMKN';
             'RTIB'; 'LTIB';
             'RANK'; 'LANK';
             'RMED'; 'LMED';
             'RHEE'; 'LHEE';
             'RTOE'; 'LTOE'};
             
                          
                          
                        
              

    [datanew] = subtract6mm(mnames,x,infile);
    mrkdata = datanew;

    done = writeTRCFile(time,mrkdata,mrknames,write_directory,infile);
    if done==1
         display(['File ' infile '_offset.trc' ' written with offset']);
    end
    
end

% end %end of 'try' for whole function

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

function done = writeTRCFile(time,mrkdata,mrknames,directory,file)
% time: nSamples vector of times
% mrkdata: nSamples x (nMarkers*3) matrix of marker xyz positions
% mrknames: nMarkers x 1 cell of marker names
% Directory you want to write to
% beginning of filename - the trc gets added automatically


if length(time)<2
    time(2)=1;
end
T=time(2)-time(1);

f=1/T;
[mk,nk]=size(mrkdata);
nk=nk/3;
fid = fopen([directory,'\',file,'.trc'],'w');
fprintf(fid,'PathFileType  4	(X/Y/Z) %s\n',directory);
fprintf(fid,'DataRate	CameraRate	NumFrames	NumMarkers	Units	OrigDataRate	OrigDataStartFrame	OrigNumFrames\n');
fprintf(fid,'%7.1f	\t%7.1f	\t%7d	\t%7d	\t mm	%7.1f	%7d	%7d\n',f,f,mk,nk,f,1,mk);
fprintf(fid,'Frame#\t');
fprintf(fid,'Time\t');
for i=1:nk
    mark=cellstr(mrknames(i,:)); mark=char(mark);
    fprintf(fid,'%s\t\t\t',mark);
end
fprintf(fid,'\n');
fprintf(fid,'		');
for i=1:nk
    if (i<10)
        fprintf(fid,'X%1d	Y%1d	Z%1d	',i,i,i);
    elseif (i<100)
        fprintf(fid,'X%2d	Y%2d	Z%2d	',i,i,i);
    else
        fprintf(fid,'X%3d	Y%3d	Z%3d	',i,i,i);
    end
end
fprintf(fid,'\n');
fprintf(fid,'\n');

for i=1:mk
    fprintf(fid,'%d',i);
    fprintf(fid,'\t%.5f',time(i));
    fprintf(fid,'\t%.3f',mrkdata(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
done=1;
end

%-------------------------------------------------------------------------%

function [datanew] = subtract6mm(mkrNames,data,infile)
%This code can be used to calculate hip joint centers based on the methods
%in Harrington et al, 2007 (J Biomech)
%Written by Cara Welker and Scott Uhlrich

% basedir = 'W:\OA_GaitRetraining\DATA\' ;
% [filename,fileDir]=uigetfile([basedir '*.trc'],'Select Trial');
% addpath(genpath('W:\OA_GaitRetraining\Matlab\common')) ;

% 
% [header data] = TRCload([fileDir filename]) ;

[timesteps,nmars] = size(data(1:end,:));
%marker indices



finalMarkerNames = {'C7'; 'CLAV';
             'RSHO'; 'LSHO';
             'RUPA'; 'LUPA';
             'RFRM'; 'LFRM';
             'RELB'; 'LELB';
             'RWRA'; 'LWRA';
             'RWRB'; 'LWRB';
             'RHJC'; 'LHJC';
             'RASI'; 'LASI';
             'RPSI'; 'LPSI';
             'RTHI'; 'LTHI';
             'RKNE'; 'LKNE';
             'RMKN'; 'LMKN';
             'RTIB'; 'LTIB';
             'RANK'; 'LANK';
             'RMED'; 'LMED';
             'RHEE'; 'LHEE';
             'RTOE'; 'LTOE'};

%marker indices
inds.C7 = (strmatch('C7',mkrNames)-1)*3+1;
C7 = data(:,inds.C7:inds.C7+2);
[nrows,ncolumns] = size(C7);
datanew = zeros(nrows,36*3);

for i = 1:36

    index = (strmatch(finalMarkerNames{i},mkrNames)-1)*3+1;
    isFull = size(index);
    if isFull == 1
        datanew(:,3*i-2:3*i) = data(:,index:index+2);
        datanew(:,3*i-1) = datanew(:,3*i-1) - 6;
    else 
        datanew(:,3*i-2:3*i) = NaN;
        display(['File ' infile '_offset.trc' ' is missing marker ' finalMarkerNames{i}]);
    end

end %for loop

end %offset and trimming markers