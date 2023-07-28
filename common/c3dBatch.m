%% Batch c3d to TRC conversion

% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% subjects = {'S0977'; 'S0978'; 'S0980'; 'S0981'; 'S0985'; 'S0989'; 'S0991'; 
%             'S0993'; 'S0995'; 'S0998'; 'S0999'; 'S1001'; 'S1002'; 'S1003';
%             'S1004'; 'S1005'; 'S1007'; 'S1009'; 'S1011'; 'S1012'; 'S1017';
%             'S1018'; 'S1019'; 'S1020'; 'S1022'; 'S1023'; 'S1024'; 'S1025';
%             'S1027'; 'S1029'; 'S1033'; 'S1034'; 'S1035'; 'S1039'; 'S1040';
%             'S1044'; 'S1047'; 'S1049'; 'S1051'; 'S1052'; 'S1053'; 'S1054'};

subjects = {'S0977'};

basedir = 'I:\Shared drives\HPL_MASPL\ProcessedData\' ;

for sub = 1:length(subjects)
    subject = subjects{sub} ;
    disp(['Analyzing Subject ' subject])

    % move to directory where this subject's files are kept
    % subjectdir = uigetdir('W:\OA_GaitRetraining\GastrocAvoidance\DATA\', 'Select the folder that contains the current subject data');
    subjectdir = [basedir subject '\'] ;

    % Go to the folder in the subject's folder where .c3d files are
    c3d_data_folder = [subjectdir '\Raw\'] ;
    names = dir(fullfile(c3d_data_folder, '*.c3d')) ;
    trialsForPrep = {names(:).name} ;
    nTrials = length(trialsForPrep);
    
    % specify where results will be printed.
    results_folder = ([subjectdir '\Raw\']);

    % Loop through the trials
    for trial= 1:nTrials

        % Get the name of the file for this trial
        c3dFile = trialsForPrep{trial};

        % Create name of trial from .c3d file name
        name = c3dFile(1:end-4) ;
        fullpath = [c3d_data_folder c3dFile] ;

        % Construct an opensimC3D object with input c3d path
        % Constructor takes full path to c3d file and an integer for forceplate
        % representation (1 = COP). 
        c3d = osimC3D(fullpath,1);

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
        if contains(name,'STATIC')
            c3d.rotateData('x',90);
            c3d.rotateData('y',0);
        else
            c3d.rotateData('x',90);
            c3d.rotateData('y',90);
        end

        % Get the c3d in different forms
        % Get OpenSim tables
        markerTable = c3d.getTable_markers();
        forceTable = c3d.getTable_forces();
        % Get as Matlab Structures
        [markerStruct, forceStruct] = c3d.getAsStructs();

        % Convert COP (mm to m) and Moments (Nmm to Nm)
        c3d.convertMillimeters2Meters();

        c3d.writeTRC([results_folder name '.trc']);
        c3d.writeMOT([results_folder name '.mot']);

        fprintf([name ' exported from .c3d to .trc & .mot' '\n']);

        % prep structure for HJC and Offset functions
        TN.basepath = results_folder;
        TN.trials = {name};

        % compute HJCs
        findHJC_RegressionOnly(TN,results_folder);

        %trim unnecessary markers and subtract 6mm from all markers 
        offsetMarkers6mm(TN,[subjectdir, 'ProcessedMARKERS']);

    end
end

