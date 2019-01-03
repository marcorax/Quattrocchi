%% m_surface
% Function decoding and interpolating disparity out of magno neuron
% activity, the interpolation is then segmented using the result of a
% previous segmentation algorithm.

function [Surface_horizontal, Surface_vertical, ScatterPerimeter] = m_surface(...
                                    OnPop, OnPopt, OffPop, OffPopt, timeL,...
                                    timeOffs, NOrient, spkThr, borders, seg_frames_l)

Surface_horizontal = zeros(size(seg_frames_l)); % Cells of horizontal disparity interpolants
Surface_vertical = zeros(size(seg_frames_l));  % Cells of vertical disparity interpolants
ScatterPerimeter = cell(length(timeL),1);

% Building the mesh grid that will be used to deal with the final results
inputsize = size(seg_frames_l);
[xq,yq]=meshgrid(1:inputsize(2),1:inputsize(1));


parfor k = 1:length(timeL)
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
    %Putting On and Off together
    tmpPos = cat(1,tmpPosOn,tmpPosOff);
    %tmpPos = tmpPosOn; using a single population for debugging    

    [UnLCenters,~,LI] = unique(tmpPos(:,1:2),'rows');
    %getting only the most firing half of the population
    UniqueLI = unique(LI);
    result = cell(2,1);
    result{1,1} = [];
    result{2,1} = [];
    for kk=1:length(UniqueLI)
        tmpDisparities=[];
        Corollas = tmpPos(LI==UniqueLI(kk),3:4);
        [UnCorollas,~,UnCorollasInd]=unique(Corollas,'rows');
        [count,values]= hist(UnCorollasInd,unique(UnCorollasInd));
        %Removing the less spiking neurons
        values = values(count>spkThr);
        count = count(count>spkThr);
        if ~isempty(count)
            tmpDisparities=sum(bsxfun(@times,...
                double(UnCorollas(values,:)-UnLCenters(UniqueLI(kk),:)),double(count)')...
                ,1)/sum(double(count));
        end
        if ~isempty(tmpDisparities)
            result{1,1}=cat(1,result{1,1},UnLCenters(UniqueLI(kk),:));
            result{2,1}=cat(1,result{2,1},tmpDisparities);
        end
    end
    % More than 2 non collinear points are needed to infer a surface in 3D, 
    % thus skip computation if less than 3 neurons are responding
    if isempty(result{1,1})
        ScatterPerimeter{k}=[];
    elseif(length(result{1,1}(:,1))>=3 && length(unique(result{1,1}(:,1)))~=1 && length(unique(result{1,1}(:,2)))~=1)
        x = double(result{1}(:,1)); 
        y = double(result{1}(:,2));
        z_h = result{2}(:,1);
        z_v = result{2}(:,2);
        interpolant_h = scatteredInterpolant(x,y,z_h, 'linear', borders);
        interpolant_v = scatteredInterpolant(x,y,z_v, 'linear', borders);
        % If the interpolants have enough points, build the surface as a 
        % 2D matrix and remove values for the background.
        Surface_horizontal(:,:,k) = interpolant_h(xq,yq);
        Surface_vertical(:,:,k) = interpolant_v(xq,yq);
        Surface_horizontal(:,:,k) = Surface_horizontal(:,:,k).*seg_frames_l(:,:,k);
        Surface_vertical(:,:,k) = Surface_vertical(:,:,k).*seg_frames_l(:,:,k);
        ScatterPerimeter{k} = [x,y,z_h,z_v];
    else % The result is not empty but it has less than 3 points, then skip
        ScatterPerimeter{k}=[];
    end
end
end






