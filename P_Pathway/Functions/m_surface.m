function [Surface_horizontal, Surface_vertical,ScatterPerimeter] = m_surface(OnPop, OnPopt,...
    OffPop, OffPopt, timeL, timeOffs, NOrient, spkThr, borders)

Surface_horizontal = cell(length(timeL)); % Cells of horizontal disparity interpolants
Surface_vertical = cell(length(timeL));  % Cells of vertical disparity interpolants

ScatterPerimeter = cell(length(timeL));

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
    if isempty(result{1,1})
        Surface_horizontal{k}=[];
        Surface_vertical{k}=[];
        ScatterPerimeter{k}=[];
    else
        x = double(result{1}(:,1)); 
        y = double(result{1}(:,2));
        z_h = result{2}(:,1);
        z_v = result{2}(:,2);
        Surface_horizontal{k} = scatteredInterpolant(x,y,z_h,'linear', borders);
        Surface_vertical{k} = scatteredInterpolant(x,y,z_v,'linear', borders);
        ScatterPerimeter{k} = [x,y,z_h,z_v];
    end
end
end






