function [] = c3dExportFunc(c3dstruct)

% Load OpenSim libs
import org.opensim.modeling.*

%initialize index variable
i = 1;
for 1:length(c3dstruct.trials)

% Get the path to a C3D file
[filename, path] = uigetfile('S1003_STATIC.c3d');
c3dpath = fullfile(path,filename);

% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP). 
c3d = osimC3D(c3dpath,1);

% Get some stats...
% Get the number of marker trajectories
nTrajectories = c3d.getNumTrajectories();
% Get the marker data rate
rMarkers = c3d.getRate_marker();
% Get the number of forces
nForces = c3d.getNumForces();
% Get the force data rate
rForces = c3d.getRate_force();

% Get Start and end time
t0 = c3d.getStartTime();
tn = c3d.getEndTime();

% Rotate the data 
c3d.rotateData('x',-90)
c3d.rotateData('y',0)

% Get the c3d in different forms
% Get OpenSim tables
markerTable = c3d.getTable_markers();
forceTable = c3d.getTable_forces();
% Get as Matlab Structures
[markerStruct forceStruct] = c3d.getAsStructs();

% Convert COP (mm to m) and Moments (Nmm to Nm)
c3d.convertMillimeters2Meters();

% Write the marker and force data to files

% Write marker data to trc file.
% c3d.writeTRC()                       Write to dir of input c3d.
% c3d.writeTRC('Walking.trc')          Write to dir of input c3d with defined file name.
% c3d.writeTRC('C:/data/Walking.trc')  Write to defined path input path.
c3d.writeTRC('S1003_STATIC_markers.trc');

% Write force data to mot file.
% c3d.writeMOT()                       Write to dir of input c3d.
% c3d.writeMOT('Walking.mot')          Write to dir of input c3d with defined file name.
% c3d.writeMOT('C:/data/Walking.mot')  Write to defined path input path.
% 
% This function assumes point and torque data are in mm and Nmm and
% converts them to m and Nm. If your C3D is already in M and Nm,
% comment out the internal function convertMillimeters2Meters()
c3d.writeMOT('S1003_STATIC_forces.mot');

i = i+1;
end

end