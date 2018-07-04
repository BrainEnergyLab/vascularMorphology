function findVascDensity(fname)
%
%written by Kira and Orla July 2018
%
%function to find the vascular density from the analysed skeleton function
%output in imageJ. Folder needs to include an excel called
%BranchInformation, a tif called CatStack, and excel files containing
%ErasedData in the title. It will output a number in metres per millimetres
%cubed, this should likely range between 0.5-3, if the numbers are way off
%this, there is probably a mistake in the preprocessing... 

%search for all the Dmap tif files across inputted directory
findExps = findFolders(fname, '*BranchInformation.csv');

for a = 1:size(findExps,2)
    
    %find excel files with info on what has been erased
    findErased = findFolders(expDir, '*ErasedData.csv');
    
    for b = 1:size(findErased,2) %loop the excel files with erased data inf
        
        %find the local exp dir
        [expDir,~]=fileparts(findExps{a});
        
        %load the excel data into a temporary variable:
        ttt = importdata(findErased{1,b});
        
        %find the index of which data column is the area value
        clear areaInd;
        areaIndCell = strfind(ttt.textdata, 'Area');
        %it automatically puts the string find info into cells, so will 
        %loop to make index a vector
        for c = 1:size(areaIndCell,2) %loop the size of index
            if isempty(areaIndCell{c})
                %the cell is empty, i.e. doesnt contain the 'area' string
                areaInd(c,:) = 0;
            else
                %contains area string, put a 1 in this index
                areaInd(c,:) = areaIndCell{c};
            end
        end %end of looping index
        clear areaIndCell;
        %find the index of which data column is the frame value
        clear frameInd;
        frameIndCell = strfind(ttt.textdata, 'Frame');
        for c = 1:size(frameIndCell,2) %loop the size of index
            if isempty(frameIndCell{c})
                %the cell is empty, i.e. doesnt contain the 'frame' string
                frameInd(c) = 0;
            else
                %contains area string, put a 1 in this index
                frameInd(c) = frameIndCell{c};
            end
        end %end of looping index
        clear frameIndCell;
        
        %get the area data & num of frames data out - save into structures:
        for c = 1:size(ttt.data,1) %loop the data rows inside excel file
            data(b).areaVal(c) = ttt.data(c,find(areaInd));
            data(b).frameNum(c) = ttt.data(c,find(frameInd));
        end %end of looping data rows
        
    end %end of looping excel files
    
    %load the ini file as we need the pixel size to convert pixels into um
    ini_file = ini2struct(cell2mat(findFolders(expDir, '*.ini')));
    %read .ini file and extract useful variables related to pixel size
    pxsz = str2double(ini_file.x_.x0x2epixel0x2esz);  % pixel size (metres)
    pxsz_um = pxsz*1000000;                        % pixel size (um)
    zstep = str2double(ini_file.x_.z0x2espacing); %z spacing in um
    
    
    %convert the erased info into volume in um
    erasedvol=0;
    for b = 1:size(data,2) %loop the excel sheets
        for c = 1:size(data(b).areaVal,2) %loop the data columns 
            %extract the area of the blob in um2 and x by frames * zstep
            %sz in microns - to give area across depth in microns, aka vol
            data(b).erasedVol(c) = data(b).areaVal(c)* ...
                (data(b).frameNum(c)*zstep);
            %also save them out not in a structure, i.e. so we can get a
            %sum across all
            if sum(erasedvol)==0
                erasedvol=data(b).erasedVol(c);
            else
                erasedvol=[erasedvol, data(b).erasedVol(c)];
            end
        end %end of loop data columns from excel
    end %end of loop excel sheets
    %sum the erased volumes to find the total erased volume
    totalerasedvolum=sum(erasedvol);
    
    %find volume of full stack - get image info from the tif file
    findStack = findFolders(expDir, '*CatStack.tif');
    infoImage=imfinfo(findStack{1,1});
    width=infoImage(1).Width*pxsz_um; %convert into um
    height=infoImage(1).Height*pxsz_um; %convert into um
    numFrames=length(infoImage)*zstep; %convert into um
    
    %find the total volume in microns of whole stack
    totalvolum=width*height*numFrames; %in um
    %convert from um cubed to mm cubed - for total and erased
    totalvolmm=totalvolum/10^9;
    totalerasedvolmm=totalerasedvolum/10^9;
    %subtract erased from total to get final vol
    finalvolmm=totalvolmm-totalerasedvolmm;
    
    %import data of branch information
    findBranchInf = findFolders(expDir, '*BranchInformation.csv');
    %load data
    clear ttt;
    ttt=importdata(findBranchInf{1,1});
    %search for branch length column
    %find the index of which data column is the area value
    clear branchInd;
    branchIndCell = strfind(ttt.textdata, 'Branch length');
    %it automatically puts the string find info into cells, so will loop to
    %make index a vector
    for c = 1:size(branchIndCell,2) %loop the size of index
        if isempty(branchIndCell{c})
            %the cell is empty, i.e. it doesnt contain the 'branch' string
            branchInd(c,:) = 0;
        else
            %contains branch string, put a 1 in this index
            branchInd(c,:) = branchIndCell{c};
        end
    end %end of looping index
    clear branchIndCell;
    
    %extract the branch lengths from the excel sheet, these are in um
    branchLength_all_um = ttt.data(:,find(branchInd));
    %convert from um to m
    branchLength_all_m = branchLength_all_um / 1000000; %m
    %sum them to get the overall branch length in z stack
    totalbranchLength = sum(branchLength_all_m);
    
    %get the final vessel density across z stack in metres per mm cubed
    vesselDensity = totalbranchLength / totalvolmm;
    
    %save the data into the exp dir
    matfile = fullfile(expDir, 'vesselDensity.mat');
    save(matfile,'vesselDensity', 'branchLength_all_m', 'totalvolmm', ...
        'totalerasedvolmm','-v7.3');
    
    
end %end of looping experiments

end %end of function
