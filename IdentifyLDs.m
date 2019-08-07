function handles = IdentifyLDs(handles)

% Help for IdentifyLDs
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Identify lipid droplets (LDs) from probabilites maps (Ilastik outputs).
% This adapted CellProfiler module resides on the CellProfiler version 1 
% see https://cellprofiler.org/previous_releases/
% and its implementation by the laboratory of Prof. Pelkmans
% see https://github.com/pelkmanslab/CellProfilerPelkmans
% 
% INPUT
% Probabilities maps need to be loaded beforehand as png in the LoadImages module
% The module IdentifyLDs should be used as a part of a image analysis pipeline containing multiple modules
%
% OUTPUT
% Labeled segmented LDs. The segmentation of the binary (black and white) image output is exported as png in the output folder.
%
% COMMENT
% For best results, the nucleus and cell outlines should be first segmented using
% implementation of the IdentifyPrimary (for the nucleus) and 
% IdentifySecondary (for the cell outline). Cells objects should be related to the
% Nucleus using the Relate module.
% To further assigned LDs, the identified LDs objects should become Children objects of the Cells (parent object)
% using the Relate module.
%
% AUTHORS:
% Christophe Freyre
%
% (c) Klemm Lab 2016

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the objects identified by this module?
%defaultVAR01 = LDs
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the intensity image that should be used for object identification?
%infotypeVAR02 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IMAGES FROM HANDLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigImage = handles.Pipeline.(ImageName);


%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%

image = OrigImage;

image = imcomplement(image);
% Create a binary image
bw = im2bw(image);

D = -bwdist(~bw,'chessboard');
D(~bw) = -Inf;

% Watershed transformation of the image
L = watershed(D,8);

labeledImage = bwlabel(L);

% Extract the properties of the object such as Area, Solidity etc
props = regionprops(labeledImage,'Area','Solidity','Perimeter','Eccentricity','Centroid','Boundingbox');

% Filter artefacts based on the measured properties
deletedobjects = find([props.Area] > 3000 & [props.Solidity] < .8); % This set of number has been empirically determined to remove artefacts of the watershed
props2 = find([props.Area] > 3000 & [props.Solidity] < .8);
delobjs = ismember(labeledImage, props2) > 0;

jjj=0;
for jjj = 1:length(deletedobjects)
    labeledImage(labeledImage==deletedobjects(jjj)) = 0;
end
deletebackground = find([props.Area] > 100000); % This number has been empirically determined and reflect the background object. It is an artefact of the watershed
jjjj=0;
for jjjj = 1:length(deletebackground)
    labeledImage(labeledImage==deletebackground(jjjj)) = 0;
end
% Label new objects before registration (see SAVE DATA TO HANDLES STRUCTURE)
NewlabeledImage = bwlabel(labeledImage);
imFinalObjects = NewlabeledImage;


%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%
%%%%%%%%%%%%%%%%%%%%%

drawnow

WatershedImage = label2rgb(L,'jet','w');
DelObjects = label2rgb(delobjs,'jet','w','shuffle');
Segmented_LDs = label2rgb(imFinalObjects,'jet','k','shuffle');

% GUI
% CPfigure(handles,'Image',ThisModuleFigureNumber);
if CPisHeadless == false
    subplot(2,2,1); CPimagesc(OrigImage,handles);
    title(['Original image',num2str(handles.Current.SetBeingAnalyzed)]);
    freezeColors
    subplot(2,2,2); CPimagesc(WatershedImage,handles);
    title(['Original Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    freezeColors
    subplot(2,2,3); CPimagesc(DelObjects,handles);
    title(['Deleted Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hold on
    freezeColors
    subplot(2,2,4); CPimagesc(Segmented_LDs,handles);
    title(['Final objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    freezeColors
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE DATA TO HANDLES STRUCTURE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldname = ['UneditedSegmented',ObjectName];%not edited for size or edge
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['SmallRemovedSegmented',ObjectName];%for IdentifySecondary.m
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['Segmented',ObjectName];%final label image
handles.Pipeline.(fieldname) = imFinalObjects;

%%% Saves location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(imFinalObjects,'Centroid');
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];   % follow CP's convention to save 0s if no object
end
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves ObjectCount, i.e. number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(imFinalObjects(:));

end

