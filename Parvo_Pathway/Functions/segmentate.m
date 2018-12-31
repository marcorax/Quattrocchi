%%% Segmentation functio
% This function is used to implement a simple surface segmentation in the 
% P_pathway, necessary to merge cells activity do compute a fine disparity
% extimation.
% The function is based on the active contours function where the contour
% is defined by the active disparity neurons coming from the Magno Pathway

function [seg_frames_l, seg_frames_r]=segmentate(FramesL, FramesR, OnPop,...
                                            OnPopt, OffPop, OffPopt,...
                                            timeL, timeOffs, NOrient)
% Number of iteration of the filling algorithm.
nIterations = 300;
seg_frames_l = cell(length(FramesL)); % Cells that will contain the segmented left frames 
seg_frames_r = cell(length(FramesR)); % Cells that will contain the segmented right frames 

    %parfor k = 1:length(FramesL)
    for k = 1:length(FramesL)
        framel= squeeze(FramesL(:,:,k));
        framer= squeeze(FramesR(:,:,k));
        tmpPosOn = [];
        tmpPosOff = [];
        for kk=1:NOrient
            %extracting events around the two frames time within a certain time
            %offsett(timeOffs)
            tmpIndOn = ((timeL(k)-timeOffs)<OnPopt{kk})&...
                (OnPopt{kk}<(timeL(k)+timeOffs));
            tmpIndOff = ((timeL(k)-timeOffs)<OffPopt{kk})&...
                (OffPopt{kk}<(timeL(k)+timeOffs));

            tmpPosOn = cat(1,tmpPosOn,OnPop{kk}(tmpIndOn,:));
            tmpPosOff = cat(1,tmpPosOff,OffPop{kk}(tmpIndOff,:));
        end
        % Putting On and Off together
        tmpPos = cat(1,tmpPosOn,tmpPosOff);
        %tmpPos = tmpPosOn; using a single population for debugging  
        % Extracting the unique centers that will be the starting point
        % for the active contours techinique
       
        % Centers for the left frame
        [UnLCenters,~,~] = unique(tmpPos(:,1:2),'rows');
        
        % Centers for the right frame
        [UnRCenters,~,~] = unique(tmpPos(:,3:4),'rows');
        
        maskL = zeros(size(framel));
        maskR = zeros(size(framer));
        
        maskL(UnLCenters) = 1;
        maskR(UnRCenters) = 1;
        % C'è un bug nella segmentazione, controlla che i centri siano
        % giusti
        tmpRes=activecontour(FramesL, maskL, nIterations);
        tmpRes=activecontour(FramesR, maskR, nIterations);
    
    end

end