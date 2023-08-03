
% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %   
% Copyright (c) 2005-2012 Stanford University and the Authors             %
% Author(s): Edith Arnold                                                 %  
%                                                                         %   
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         % 
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
% ----------------------------------------------------------------------- % 

% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% subjects = {'S0977'; 'S0978'; 'S0980'; 'S0981'; 'S0985'; 'S0989'; 'S0991'; 
%             'S0993'; 'S0995'; 'S0998'; 'S0999'; 'S1001'; 'S1002'; 'S1003';
%             'S1004'; 'S1005'; 'S1007'; 'S1009'; 'S1011'; 'S1012'; 'S1017';
%             'S1018'; 'S1019'; 'S1020'; 'S1022'; 'S1023'; 'S1024'; 'S1025';
%             'S1027'; 'S1029'; 'S1033'; 'S1034'; 'S1035'; 'S1039'; 'S1040';
%             'S1044'; 'S1047'; 'S1049'; 'S1051'; 'S1052'; 'S1053'; 'S1054'};

% masses = [72.2; 105.9; 78.9; 65.0; 76.1; 73.5; 64.7; 
%           65.8; 58.3; 52.7; 61.9; 65.8; 39.8; 65.8;
%           59.9; 65.8; 73.1; 60.6; 65.8; 79.8; 80.6;
%           102.3; 65.8; 47.3; 46.5; 58.7; 69.2; 70.1;
%           49.2; 70.7; 58.8; 71.2; 70.5; 82.3; 105.8;
%           60.2; 49.5; 66.4; 73.7; 74.0; 46.4; 65.3];


subjects = {'S1003'};
masses = [65.8];

basedir = 'C:\SharedGDrive\HPL_MASPL\ProcessedData\' ;

% modelName = 'LaiUhlrich2022_marked.osim';
modelName = 'LaiUhlrich2022.osim';
markerSetFileName = 'C:\MyRepositories_Xander\opencap-core_ACL\opensimPipeline\Models\Vicon_markers.xml';

batchScaleSettingsFileName = 'C:\MyRepositories_Xander\opencap-core_ACL\opensimPipeline\Scaling\Setup_scaling_Vicon.xml' ;

for sub = 1:length(subjects)
    subject = subjects{sub} ;
    disp(['Analyzing Subject ' subject])

    % move to directory where this subject's files are kept
    % subjectdir = uigetdir('W:\OA_GaitRetraining\GastrocAvoidance\DATA\', 'Select the folder that contains the current subject data');
    subjectdir = [basedir subject '\'] ;

    % Go to the folder in the subject's folder where .trc files are
    trc_data_folder = [subjectdir 'ProcessedMARKERS\'] ;
    names = dir(fullfile(trc_data_folder, '*STATIC.trc')) ;
    trialsForScale = {names(:).name} ;
    nTrials = length(trialsForScale);
    
    % specify where results will be printed.
    results_folder = ([subjectdir]);

    % Get and operate on the files
    % Choose a generic setup file to work from
    scaleTool = ScaleTool(batchScaleSettingsFileName);
    modelScaler = scaleTool.getModelScaler;
    markerPlacer = scaleTool.getMarkerPlacer;
    genericModelMaker = scaleTool.getGenericModelMaker;
    % Get the model
    % Load the model and initialize
    model = Model([basedir modelName]);
    model.initSystem();

    % Tell Tool to use the loaded model
    % scaleTool.setModel(model);

    % Loop through the trials
    for trial= 1:nTrials;

        %Generic setup file to work from
        %ikTool = InverseKinematicsTool(['S:\Human Performance Lab\MotionAnalysis\OpenSim Models and Files\Generic_ACL_IKSetup_report_markers.xml']);
        %model = Model([modelFilePath modelFile]);
        %model.initSystem();
        %ikTool.setModel(model);

        % Get the name of the file for this trial
        markerFile = trialsForScale{trial};

        % Create name of trial from .trc file name
        name = markerFile(1:end-4) ;
        fullpath = [trc_data_folder markerFile] ;

        % Get trc data to determine time range
        markerData = MarkerData(fullpath);

        % Get initial and final time 
        initial_time = markerData.getStartFrameTime();
        final_time = markerData.getLastFrameTime();
        timerange = ArrayDouble(0, 2);
        timerange.setitem(0,initial_time);
        timerange.setitem(1,final_time);


        % Setup the scaleTool for this trial
        scaleTool.setName([modelName(1:end-5) '_' subject '.osim']);
        modelScaler.setMarkerFileName(fullpath);
        markerPlacer.setMarkerFileName(fullpath);
        modelScaler.setTimeRange(timerange);
        markerPlacer.setTimeRange(timerange);
        scaleTool.setSubjectMass(masses(sub));
        genericModelMaker.setModelFileName([basedir modelName]);
        genericModelMaker.setMarkerSetFileName(markerSetFileName);
        modelScaler.setOutputModelFileName([results_folder modelName(1:end-5) '_' subject '.osim']);
        markerPlacer.setOutputModelFileName([results_folder modelName(1:end-5) '_' subject '.osim']);
        markerPlacer.setOutputMarkerFileName([results_folder subject '_markers.xml']);

        % scaleTool.createModel();
        % scaleTool.createModel([results_folder 'LaiUhlrich2022_marked_scaled_' subject '.osim']);

        % Save the settings in a setup file
        outfile = ['Setup_Scale_' name '.xml'];
        try
            scaleTool.print([results_folder outfile]);
        catch
            mkdir([results_folder name])
            scaleTool.print([results_folder outfile]);
        end
        % warning off
        % mkdir([results_folder name '\output\'])
        % warning on

        fprintf(['Performing Scaling for ' name '\n']);
        % Run ScaleTool in the command line
        % scaleTool = ScaleTool([results_folder outfile]); 
        % scaleTool.run()
        myCommand = ['opensim-cmd run-tool ' results_folder outfile];
        system(myCommand)

        %Give unique name to marker location .sto file (if set to 'true' in
        %setup file)
%         movefile([name '_ik_model_marker_locations.sto'],[results_folder name '\output\ik_model_marker_locations.sto'])
%         movefile([name '_ik_marker_errors.sto'],[results_folder name '\output\ik_model_marker_errors.sto'])

        %clear ikTool;

    end
end