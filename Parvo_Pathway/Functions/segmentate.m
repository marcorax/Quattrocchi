%%% Segmentation functio
% This function is used to implement a simple surface segmentation in the 
% P_pathway, necessary to merge cells activity do compute a fine disparity
% extimation.
% The function is based on the active contours function where the contour
% is defined by the active disparity neurons coming from the Magno Pathway

function [seg_frames]=segmentate(Frames, OnPop, OnPopt, OffPop, OffPopt,...
                                                timeL, timeOffs, NOrient,...
                                                nIterations, bord_thick)

seg_frames = zeros(size(Frames)); % Cells that will contain the segmented left frames 

    parfor k = 1:length(Frames)
        frame= squeeze(Frames(:,:,k));
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
        [UnCenters,~,~] = unique(tmpPos(:,1:2),'rows');
        
        mask = zeros(size(frame));
        
        % In matlab I need to use 1D indeces to access a predefined list of
        % point in a matrix
        Index=sub2ind(size(frame),UnCenters(:,2),UnCenters(:,1));
        
        % Maximum brightness value for a 16 bit unsigned value (The same
        % of the sensor used: DAVIS240)
        mask(Index) = 1;
        
        if(bord_thick~=1)
            % Increase the size of the mask borders
            mask = logical(conv2(mask, ones(bord_thick), 'same'));
        end
        
        tmp_seg_frame=double(activecontour(frame, mask, nIterations));
        % putting zeros to NaN
        tmp_seg_frame(tmp_seg_frame == 0)= NaN;
        seg_frames(:,:,k) = tmp_seg_frame;
    end

end