
function findDiamDistMap(fname)
%
%function written by Kira & Orla, June 2018
%
%function loads in the diameter map tif file (from ImageJ) and the dmap
%coordinates (from excel, exported from imagej) to work out the vessel
%diameter across each skeleton/branch from a zstack 
%the diameters outputted are in um, so the user needs to make sure that the
%image has been SCALED in imageJ for these to be accurate
%
%INPUTS-
%fname = top dir which contains the 3Dcoords.xlsx files, note each dir with 
%an excel must also have a tif file called Dmap.tif 
%
%OUTPUTS-
%nothing outputted from function, but will save a .mat file called
%extractedDiameters - with the data processed in cohesive matlab form (i.e.
%saved in separate cells for each skeleton and branch)
%
%functions to run before this:
%NEED to run get3DSkelCoords.bsh (in imageJ) before and save 3Dcoords.xlsx
%
%other functions needed in the path:
%findFolders

%search for all the Dmap tif files across inputted directory
findExps = findFolders(fname, '*3Dcoords.xlsx');

%loop the exp folders which contain the relevant files
for a = 1:size(findExps,2)

%find the local exp dir
[expDir,dmapFile]=fileparts(findExps{a});

%% load the tif file with the processed images:
%find the dmap file which to get the image dimension info
%NB/ will not be using pixel size, as images should already be scaled
chinfo = cell2mat(findFolders(expDir, '*Dmap.tif'));
%get image info:
info = imfinfo(chinfo);
width = info(1).Width;    %width of frame (pixels)
height = info(1).Height;  %height of frame (pixels)
numFrames = numel(info); %number of frames (frames)

disp('loading DMap tif file...'); 
%load the dmap image tif file into the workspace
%allocate zeros for Dmap import. Will need to edit if not 256 x 256
DMap = zeros(width, height, numFrames);
%import DMap tiff into empty matrix
for k = 1:numFrames
    %Distance map TIFF
     DMap(:,:,k) = imread(findExps{1,a}, k);
end

%% load the 3D coordinates and put into organised cells:
%load the excel file with the 3D skel coords in - as imported from imageJ
rawCoordData = importdata('3Dcoords.xlsx');

%find all the skeleton labels, so can keep vars separated
findSkelStr = strfind(rawCoordData,'Skeleton');
for b=1:size(findSkelStr,1)
    if isempty(findSkelStr{b,1})
        findSkelStr{b,1}=0;
    end
end
skelInd=find(cell2mat(findSkelStr)>=1);
%find all the branch labels, so can keep vars separated by branch too 
findBranchStr = strfind(rawCoordData,'Branch');
for b=1:size(findBranchStr,1)
    if isempty(findBranchStr{b,1})
        findBranchStr{b,1}=0;
    end
end
branchInd=find(cell2mat(findBranchStr)>=1);

%give the branchInd a skelInd label - for looping branches within skel
for b = 1:size(skelInd,1) %loop skel ind
    if b~= size(skelInd,1)
        for c = 1:size(branchInd,1) %loop branchInd
            %check which skelInd switch point it falls into
            if branchInd(c) > skelInd(b) && branchInd(c) < skelInd(b+1)
                %give it the skel label
                branchInd(c,2) = b;
            end
        end %end of check if last skel pt
    else
        %if it is the last skel pt, give it the last label
        branchInd(c,2) = size(skelInd,1);
    end %end of loop branch ind
end %end of loop skel ind

%take all the branches inside the skeletons, and put in separate cells - so
%data is organised in a more cohesive manner
disp('sorting data by skeletons and branches...'); %inform user
for b = 1:size(skelInd,1) %loop skeleton Ind
    
    %find how many branches within this skel for the num of cells
    cellSzInd=find(branchInd(:,2)==b);
    %so know how many cells for 2nd dim
    for c = 1:size(cellSzInd,1) %loop branches
        
            %find the index for the start and end of the branch
            start=branchInd(cellSzInd(c),1);
            %calculate stop pt from index, unless last index switch pt,
            %then need to use the size of the data
            if b~=size(skelInd,1)%check if last ind pt
                %calc stop pt from next ind pt - 1
                stop=branchInd(cellSzInd(c)+1,1)-1;
            else
                %if last ind pt, use size of data
                stop=size(rawCoordData,1);
            end %end of check if last skel switch pt
            %loop the data pts within the branch and put into temp matrix
            clear ttt;
            for d = 1:stop-start
                if d==1 %if it is the first data pt within the branch
                    ttt=str2num(cell2mat(rawCoordData(start+d)));
                else %later data in branch
                    ttt=[ttt; str2num(cell2mat(rawCoordData(start+d)))];
                end %check which data pt in branch
            end %end of loop data pts
        
        %put data into correct skel pt and branch pt cell
        rawCoordDataSorted{b}{c}=ttt;
        
    end %end of loop branches
end %end of loop skels

%% extract the dmap value and diam for each xyz coord

%loop through the coordinates and check the image
for b = 1:size(rawCoordDataSorted,2) %loop the skeletons
    for c = 1:size(rawCoordDataSorted{b},2) %loop the branches
        
        %extract the raw coordinate values from index created, so can find
        %the dmap value from image
        clear xyzcoord_ttt;
        xyzcoord_ttt=rawCoordDataSorted{b}{c};
        
        for d = 1:size(xyzcoord_ttt,1)
            %x and y are flipped in matlab images, because they go col then
            %row, rather than row then col as in imagej
            %so ask for coord 2 then coord 1
            DMapCoordSortedVal{b}{c}(d,:) = ...
                DMap(xyzcoord_ttt(d,2),xyzcoord_ttt(d,1),xyzcoord_ttt(d,3));
            %NB you can use this number to check if it matches the value in
            %matlab for the corresponding coords
            
            %get diameter by multiplying by 2
            DMapCoordSortedDiam{b}{c}(d,:) = 2* ... 
                DMap(xyzcoord_ttt(d,2),xyzcoord_ttt(d,1), ...
                xyzcoord_ttt(d,3));
            
        end
        
    end %end of branches loop
end %end of skel loop

%% get average diameters 

%get average for each branch
for b = 1:size(DMapCoordSortedDiam,2) %loop skeletons
    for c = 1:size(DMapCoordSortedDiam{b},2)%loop branches
        %average across branch
        DMapBranchDiam{b}{c}=nanmean(DMapCoordSortedDiam{b}{c});
    end %end of loop branches
end %end of loop skeletons

%% save data into mat file in exp dir
matfile = fullfile(expDir, 'extractedDiameters.mat');
save(matfile,'DMapBranchDiam','DMapCoordSortedDiam','DMap', ...
    'DMapCoordSortedVal','rawCoordData','rawCoordDataSorted', '-v7.3');


end %end of looping exp dirs

end %end of function 
