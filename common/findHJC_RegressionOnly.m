function [] = findHJC_RegressionOnly(trialNames,outputPath)
% Created by Xander Sahin
% This program is modified from the findHJC_Dynamic_Regression function
% created by Amy Silder, Julie Kolesar, and Scott Ulrich. It computes 
% hip joint centers based on regression only. New marker positions are 
% written to the *.trc file with RHJC & LHJC. 
% 
% Inputs:
% trialNames is a structure with the following fields:
%   basepath: path to where *.trc files are
%   trials: a cell of relative paths to all the trials to which HJCs are to
%   be added
%     
% outputPath: fullpath of folder to which new files get written

%
% Output:
% None directly, new TRC files get written

% try %whole main function

%file and folder naming    
% userDir = getuserdir ;
% tempOutputPath = [userDir '\Temp\Sub' num2str(round(rand*1000000)) '\'] ;
finalFolderName = '\TRCs_w_HJCs\';
outputPathFull = [outputPath,finalFolderName] 

%Select files to which the HJC locations should be added (multiple must be
%selected).  The same directory will be used to create a file for new files
%with HJC locations.
files = trialNames.trials ;
directory = trialNames.basepath ;

files=strvcat(files(:,:));

if length(files)==1;
    display('No Files Selected to Write HJC.  Try Again.');
    return;
end

directory=strvcat(directory);
infile=files(1,:); inpath=directory;


%Reformat file list and directory information

[a b]=size(files);
write_directory=outputPathFull;

%read each *.trc file and append HJC markers
for j=1:a;
    infile=files(j,:);
    inpath=directory;
    x=0; mnames=0; tx=0; marks=0; time=0;
    [x,tx,sfx,nsx,nmrk,mnames,file,inpath]=load_trc(['\' infile],inpath);   
    
    %compute time from the sampling frequency
    [a,b]=size(x);
    for i=1:a;
        time(i,1)=i/sfx-1/sfx;
    end
    nfile=strcat(files(j,:));
    filen=nfile(1,1:(length(nfile)-4));
     
    % Add regression based HJCs
    mrknames=char([mnames;cellstr('RHJC');cellstr('LHJC')]);     
    [datanew] = regressionHJC(mnames,x);
    mrkdata = [x,datanew];

    done = writeTRCFile(time,mrkdata,mrknames,write_directory,filen);
    if done==1
         display(['File ' infile ' written with HJC locations']);
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
if (n==0)
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
elseif (n==1)
    file = infile(1:length(infile)-4);
    fid=fopen(infile,'r');
else (n==1)
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
pos=pos/1000;

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
% fid = fopen([directory,'/',file,'.trc'],'w');
fid = fopen([directory,'/',file],'w');
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

function [datanew] = regressionHJC(mkrNames,data)
%This code can be used to calculate hip joint centers based on the methods
%in Harrington et al, 2007 (J Biomech)
%Written by Cara Welker and Scott Uhlrich

% basedir = 'W:\OA_GaitRetraining\DATA\' ;
% [filename,fileDir]=uigetfile([basedir '*.trc'],'Select Trial');
% addpath(genpath('W:\OA_GaitRetraining\Matlab\common')) ;

pelvisMarkNames = {'RASI','RPSI','LASI','LPSI'} ;
% 
% [header data] = TRCload([fileDir filename]) ;

[timesteps,nmars] = size(data(1:end,:));
%marker indices
inds.R_ASIS = (strmatch(pelvisMarkNames{1},mkrNames)-1)*3+1 ;
inds.R_PSIS = (strmatch(pelvisMarkNames{2},mkrNames)-1)*3+1 ;
inds.L_ASIS = (strmatch(pelvisMarkNames{3},mkrNames)-1)*3+1 ;
inds.L_PSIS = (strmatch(pelvisMarkNames{4},mkrNames)-1)*3+1 ;

L_PSIS = data(:,inds.L_PSIS:inds.L_PSIS+2)*1000;
R_PSIS = data(:,inds.R_PSIS:inds.R_PSIS+2)*1000;
L_ASIS = data(:,inds.L_ASIS:inds.L_ASIS+2)*1000;
R_ASIS = data(:,inds.R_ASIS:inds.R_ASIS+2)*1000;

[nrows,ncolumns] = size(L_ASIS);

%new HJC values
L_HJC = zeros(nrows,3);
R_HJC = zeros(nrows,3);

%origin halfway between two asis markers
for ii = 1:nrows
    PW = norm(L_ASIS(ii,:) - R_ASIS(ii,:));         %pelvic width
    
    mid_PSIS = (L_PSIS(ii,:) + R_PSIS(ii,:))./2;
    mid_ASIS = (L_ASIS(ii,:) + R_ASIS(ii,:))./2;
    origin = mid_ASIS;
    PD = norm(mid_ASIS - mid_PSIS);      %pelvic depth
    
    z = (R_PSIS(ii,:) - L_PSIS(ii,:))/norm(R_PSIS(ii,:) - L_PSIS(ii,:));
    plane = cross(z,(R_ASIS(ii,:)-mid_PSIS));
    plane_norm = plane/norm(plane);
    x = cross(plane_norm,z);
    y = cross(z,x);
    
    
    R_HJC(ii,1:3) = origin+(-0.24*PD-9.9)*x  ...
        + (-0.3*PW-10.9)*y ...
        + (0.33*PW+7.3)*z;
    
    L_HJC(ii,1:3) = origin+(-0.24*PD-9.9)*x ... 
        + (-0.3*PW-10.9)*y ...
        + (-0.33*PW+7.3)*z;
end
datanew = [R_HJC L_HJC]/1000;
% headernew = header;
% headernew.markername = horzcat(header.markername,...
%     {'RHJC_reg','','','LHJC_reg','',''});
% 
% time = data(:,strmatch('Time',header.markername));
% mrkdata = datanew(:,3:end) ; 
% mrknames_new = headernew.markername(3:end)' ;
% mrknames_new2 = cell(length(mrknames_new)/3,1) ;
% 
% j = 1 ;
% for i = 1:length(mrknames_new)
%     cellVal = mrknames_new{i} 
%     ~strcmp(cellVal,'')
%     if ~strcmp(cellVal,'')
%     mrknames_new2{j}  = cellVal ;
%     j = j+1 ;
%     end
% end
end %regression HJC