%%% Segmentation functio
% This function is used to implement a simple surface segmentation in the 
% P_pathway, necessary to merge cells activity do compute a fine disparity
% extimation.
% The function is based on the active contours function where the contour
% is defined by the active disparity neurons coming from the Magno Pathway

function [seg_frames_l,seg_frames_r]=segmentate(FramesL, FramesR, OnPop,...
                                                OnPopt, OffPop, OffPopt,...
                                                timeL, timeOffs, NOrient,...
                                                nIterations, bord_thick)

seg_frames_l = zeros(size(FramesL)); % Cells that will contain the segmented left frames 
seg_frames_r = zeros(size(FramesL)); % Cells that will contain the segmented right frames 
% If required disable warnings (likely related to activecontour)


    parfor k = 1:length(FramesL)
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
       
        % Unique centers for the left and rigth frames
        [UnCentersL,~,~] = unique(tmpPos(:,1:2),'rows');
        [UnCentersR,~,~] = unique(tmpPos(:,3:4),'rows');
        maskL = zeros(size(framel));
        maskR = zeros(size(framer));

        % In Matlab I need to use 1D indeces to access a predefined list of
        % point in a matrix
        IndexL=sub2ind(size(framel),UnCentersL(:,2),UnCentersL(:,1));
        IndexR=sub2ind(size(framer),UnCentersR(:,2),UnCentersR(:,1));

        maskL(IndexL) = 1;
        maskR(IndexR) = 1;

        if(bord_thick~=1)
            % Increase the size of the mask borders
            maskL = logical(conv2(maskL, ones(bord_thick), 'same'));
            maskR = logical(conv2(maskR, ones(bord_thick), 'same'));
       end
        
        tmp_seg_frame_l=double(activecontour(framel, maskL, nIterations));
        tmp_seg_frame_r=double(activecontour(framer, maskR, nIterations));
       
        % putting zeros to NaN
        tmp_seg_frame_l(tmp_seg_frame_l == 0)= NaN;
        tmp_seg_frame_r(tmp_seg_frame_r == 0)= NaN;

        seg_frames_l(:,:,k) = tmp_seg_frame_l;
        seg_frames_r(:,:,k) = tmp_seg_frame_r;

    end

end