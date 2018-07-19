
PlotResult = 0;
FullLat = [-90:gstep:90];
FullLon = 0:gstep:360-gstep;

% compile daily stratospheric columns
if Mn < 10
    if Windows == 0
        Files = dir(sprintf('%sMFiles/%s/no2track%d0%d*.mat',DailyStratDir,Instrument,Year,Mn));
    else
        Files = dir(sprintf('%sMFiles\\%s\\no2track%d0%d*.mat',DailyStratDir,Instrument,Year,Mn));        
    end
else
    if Windows == 0
        Files = dir(sprintf('%sMFiles/%s/no2track%d%d*.mat',DailyStratDir,Instrument,Year,Mn));
    else
        Files = dir(sprintf('%sMFiles\\%s\\no2track%d%d*.mat',DailyStratDir,Instrument,Year,Mn));        
    end
end

if ~isempty(Files)
    
    strat_2d = zeros(length(FullLat),length(FullLon));
    strat_2d_area = zeros(length(FullLat),length(FullLon));
    
    for i = 1:length(Files)
        if Windows == 0
            load(sprintf('%sMFiles/%s/%s',DailyStratDir,Instrument,Files(i).name))
        else
            load(sprintf('%sMFiles\\%s\\%s',DailyStratDir,Instrument,Files(i).name))            
        end
        
        % shift from bottom left to centre (obsolete (April 2008))
        %latgrid = latgrid + gstep/2;
        %longrid = longrid + gstep/2;
        
        longrid(longrid<-1e-6) = longrid(longrid<-1e-6) + 360;
        scdstr(area==0) = 0;
        
        % place in overall grid
        a = find(abs(latgrid(1)-FullLat)<1e-6);
        b = find(abs(longrid(1)-FullLon)<1e-6);
        
        zb = find(abs(longrid)<1e-6);
        if ~isempty(zb) & zb~=1
            strat_2d(a:a+length(latgrid)-1,b:end) = strat_2d(a:a+length(latgrid)-1,b:end) + scdstr(:,1:zb-1).*area(:,1:zb-1);
            strat_2d_area(a:a+length(latgrid)-1,b:end) = strat_2d_area(a:a+length(latgrid)-1,b:end) + area(:,1:zb-1);

            strat_2d(a:a+length(latgrid)-1,1:length(longrid)-zb+1) = strat_2d(a:a+length(latgrid)-1,1:length(longrid)-zb+1) + scdstr(:,zb:end).*area(:,zb:end);
            strat_2d_area(a:a+length(latgrid)-1,1:length(longrid)-zb+1) = strat_2d_area(a:a+length(latgrid)-1,1:length(longrid)-zb+1) + area(:,zb:end);
            
        else
            strat_2d(a:a+length(latgrid)-1,b:b+length(longrid)-1) = strat_2d(a:a+length(latgrid)-1,b:b+length(longrid)-1) + scdstr.*area;
            strat_2d_area(a:a+length(latgrid)-1,b:b+length(longrid)-1) = strat_2d_area(a:a+length(latgrid)-1,b:b+length(longrid)-1) + area;
        end
        disp(sprintf('get_strat_2d: %s read in',Files(i).name));
    end
    
    strat_2d = strat_2d./strat_2d_area;
    
    strat_1d = get_strat_1d(FullLat,FullLon,strat_2d);
    
    for i = 1:size(strat_2d,1)
        strat_scale_2d(i,:) = strat_2d(i,:)./strat_1d(i);
    end
else
    
    disp('No 2D Stratospheric Files! Using 1D')
    strat_scale_2d = ones(length(FullLat),length(FullLon));
    
end

% crop scalar to original size
strat_scale_2d = strat_scale_2d(find(Lat(1)==FullLat):find(Lat(end)==FullLat),:);

s = strat_scale_2d;
NanFind = isnan(s);
f = [5 10];
while(~sum(sum(NanFind))==0)
    [spot spot2] = find(NanFind==1);
    for i = 1:length(spot)
        x = max(spot(i)-f(1),1);
        x2 = min(spot(i)+f(1),size(s,1));
        y = max(spot2(i)-f(2),1);
        y2 = min(spot2(i)+f(2),size(s,2));
        
        s(spot(i),spot2(i)) = nanmedian(nanmedian(s(x:x2,y:y2)));
    end
    NanFind = isnan(s);
end
h = fspecial('disk',10);
s = imfilter(s,h,'circular');

strat_scale_2d = s;


if PlotResult
	plotlat = [];
	for i = 1:length(Lon)-1
        plotlat = [plotlat (Lat-gstep/2)'];
	end
	
	plotlon = [];
	for i = 1:length(Lat)
        plotlon = [plotlon; (Lon(1:end-1)-gstep/2)];
	end
	
	figure
    
    h(1) = worldmap('world');
    setm(gca,'flatlimit',[-50 60]);
	surfm(plotlat,plotlon,strat_scale_2d)
    Kids = get(h(1),'Children');
    for i = 1:size(Kids,1)
        if strcmp(get(Kids(i),'Tag'),'PLabel') | strcmp(get(Kids(i),'Tag'),'MLabel')
            set(Kids(i),'Visible','off');
        end
    end       
	set(gca,'clim',[0.8 1.2]);
    
	h(2) = colorbar('horiz');
    h(3) = title('Stratospheric Scaler');
    set(h(1),'Position',[0.01 0.25 0.98 0.67], 'visible','off')
	set(h(2),'Position',[0.13 0.32 0.7750 0.03])

    s = strat_scale_2d;
    NanFind = isnan(s);
    f = [5 10];
    while(~sum(sum(NanFind))==0)
        [spot spot2] = find(NanFind==1);
        for i = 1:length(spot)
            x = max(spot(i)-f(1),1);
            x2 = min(spot(i)+f(1),size(s,1));
            y = max(spot2(i)-f(2),1);
            y2 = min(spot2(i)+f(2),size(s,2));
            
            s(spot(i),spot2(i)) = nanmedian(nanmedian(s(x:x2,y:y2)));
        end
        NanFind = isnan(s);
    end
    h = fspecial('disk',10);
    surfm(plotlat,plotlon,imfilter(s,h,'circular'))
    %surfm(plotlat,plotlon,s)
end
