%Get monthly mean tropospheric vertical column

if sum(gstep == [2 2.5]) == 2
	lat = [-51:gstep(1):67];
	lon = -1.25:gstep(2):358.75;
    nexus = [lat(4) lon(4)];
elseif sum(gstep == [1 1]) == 2
	lat = [-50:gstep(1):66];
	lon = 0:gstep(2):360-gstep(2); 
    nexus = [lat(4) lon(4)];    
elseif sum(gstep == [1 1.25]) == 2
	lat = [-50:1:70];
	lon = 0:gstep(2):360-gstep(2); 
    nexus = [lat(4) lon(4)];    
else
	lat = [-50:gstep(1):66];
	lon = 0:gstep(2):360-gstep(2); 
    nexus = [lat(4) lon(4)]; % [gstep(1)/2 gstep(2)/2];
end

% Retrieve list of files
if strcmp(Instrument,'GOME')
    if Windows == 1
        Files = dir(sprintf('%s%d-%s%s\\*%s',OutDir,Year,MonthName(Mn,:),DirExt,FileExt));
    else
        Files = dir(sprintf('%s%d-%s%s/*%s',OutDir,Year,MonthName(Mn,:),DirExt,FileExt));
    end
elseif strcmp(Instrument,'SCIAMACHY')
    if strcmp(Source,'Dal')
        if amfv58 == 1
            Files = dir(sprintf('%s%d-%s%s/AMF-orb*',OutDir,Year,MonthName(Mn,:),DirExt));            
        else
            if Windows == 1
                Files = dir(sprintf('%s%d-%s%s\\orb*',OutDir,Year,MonthName(Mn,:),DirExt));
            else
                Files = dir(sprintf('%s%d-%s%s/orb*',OutDir,Year,MonthName(Mn,:),DirExt));
            end
        end
    elseif strcmp(Source,'KNMI')
        if Windows == 1
            Files = dir(sprintf('%s%d-%s%s\\SCIA-*',OutDir,Year,MonthName(Mn,:),DirExt));
        else
            Files = dir(sprintf('%s%d-%s%s/SCIA-*',OutDir,Year,MonthName(Mn,:),DirExt));
        end        
    end
    weekcnt = zeros(length(lat),length(lon));
    weekmno2 = zeros(length(lat),length(lon));
	weeksno2 = zeros(length(lat),length(lon));
end

cnt = zeros(length(lat),length(lon));
amf = zeros(length(lat),length(lon));
mno2 = zeros(length(lat),length(lon));
sno2 = zeros(length(lat),length(lon));
wcld = zeros(length(lat),length(lon));
err = zeros(length(lat),length(lon));
sc = zeros(length(lat),length(lon));
stratdata = zeros(length(lat),length(lon));

for i = 1:size(Files,1)

    if Windows == 1
        dat = load(sprintf('%s%d-%s%s\\%s',OutDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
    else
        dat = load(sprintf('%s%d-%s%s/%s',OutDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
    end
    
    if (length(dat)>0)
        if ExcludeNegatives == 1
            dat(dat(:,17) < 0,:) = [];
            ExNeg = '-ExNeg';
        else
            ExNeg = '';
        end
    end
    
    if (length(dat)>0)
        if Frescov5 == 0
            dat = dat((dat(:,13)<65 & dat(:,13)>-49),:); %filter for latitude
            dat = dat((dat(:,4)<0.5),:); %filter for fraction of backscatterd refl
        else
            % no format difference in outputted files
            dat = dat((dat(:,13)<65 & dat(:,13)>-49),:); %filter for latitude
            dat = dat((dat(:,4)<0.5),:); %filter for fraction of backscatterd refl            
        end
    end
    
    if length(dat)>0

        %[ave,area,latgrid,longrid]=realign(dat,gstep(1),gstep(2),[1],[2 6; 3 7; 4 8; 5 9],nexus);
        %[ave,area,latgrid,longrid]=realignv2(dat,gstep(1),gstep(2),[15,16,17,4,3,2,20],[5 6; 7 8; 9 10; 11 12],nexus);
        [ave,area,latgrid,longrid]=realignv2(dat,gstep(1),gstep(2),[15,16,17,4,3,2],[5 6; 7 8; 9 10; 11 12],nexus);
        
        %disp(sprintf('%d: min:%.2f, max: %.2f',i,min(longrid),max(longrid)));
        
        ave(isnan(ave)) = 0;
        longrid(longrid>max(lon)) = longrid(longrid>max(lon)) - 360;
        longrid(longrid<min(lon)) = longrid(longrid<min(lon)) + 360;
        
        % find indices of new data
        [a ind] = min(abs(lat-latgrid(1)));
        ilat = ind:ind+length(latgrid)-1;
        
		%if length(longrid)>length(lon)
        %    disp(sprintf('in the mysterious condition i = %d',i));
        %    beep; pause(0.5); beep; pause(0.5); beep;
	    %		area=area(:,1:length(longrid));
		% ave=ave(:,1:length(longrid),:);
		%end

        % find indices of new data
        %[a ind] = min(abs(lat-latgrid(1)));
        %ilat = ind:ind+length(latgrid)-1;
        
        %[a ind] = min(abs(lon-longrid(1)));
		%if (longrid(length(longrid))) > longrid(1) 
    	%	ilon = ind:ind+length(longrid)-1;
		%else
            %disp(sprintf('data is wrapped! i = %d',i));
            %beep; pause(0.5); beep; pause(0.5); beep;
    	%	ilon = [ind:length(lon) 1:length(longrid)-(length(lon)-ind)-1];
		%end %if wrap

        %if max(ilon) <= size(cnt,2)
        
        for loncycle = 1:size(area,2)
           
            [a ilon] = min(abs(lon - longrid(loncycle)));
        
            cnt(ilat,ilon)=cnt(ilat,ilon)+area(:,loncycle);
            amf(ilat,ilon)=amf(ilat,ilon)+ave(:,loncycle,1).*area(:,loncycle);
            mno2(ilat,ilon)=mno2(ilat,ilon)+ave(:,loncycle,2).*area(:,loncycle);%+ave(:,:,2).*area;
            sno2(ilat,ilon)=sno2(ilat,ilon)+ave(:,loncycle,3).*area(:,loncycle);%+ave(:,:,3).*area;
            wcld(ilat,ilon)=wcld(ilat,ilon)+ave(:,loncycle,4).*area(:,loncycle);%+ave(:,:,4).*area;
            err(ilat,ilon)=err(ilat,ilon)+ave(:,loncycle,5).*area(:,loncycle);%+ave(:,:,5).*area;
            sc(ilat,ilon)=sc(ilat,ilon)+ave(:,loncycle,6).*area(:,loncycle);%+ave(:,:,6).*area;
            if size(ave,3) > 6
                stratdata(ilat,ilon)=stratdata(ilat,ilon)+ave(:,loncycle,7).*area(:,loncycle);%+ave(:,:,6).*area;
            end

            if strcmp(Instrument, 'SCIAMACHY')
                if strcmp(Source,'Dal')
                    if amfv58 == 1
                        day = weekday(str2num(Files(i).name(19:26)));
                    else
                        day = weekday(str2num(Files(i).name(15:22)));
                    end
                elseif strcmp(Source,'KNMI')
                    day = weekday(str2num(Files(i).name(6:9))*10000 + MonthNameToNum(Files(i).name(11:13))*100 + str2num(Files(i).name(15:16)));                
                end
                if ~(strcmp(day,'Saturday') | strcmp(day,'Sunday'))
                    %weekcnt(ilat,ilon)=weekcnt(ilat,ilon)+area;
                    %weekmno2(ilat,ilon)=weekmno2(ilat,ilon)+ave(:,:,2).*area;
                    %weeksno2(ilat,ilon)=weeksno2(ilat,ilon)+ave(:,:,3).*area;
                    weekcnt(ilat,ilon)=weekcnt(ilat,ilon)+area(:,loncycle);
                    weekmno2(ilat,ilon)=weekmno2(ilat,ilon)+ave(:,loncycle,2).*area(:,loncycle);
                    weeksno2(ilat,ilon)=weeksno2(ilat,ilon)+ave(:,loncycle,3).*area(:,loncycle);
                end
            end
        end
    end
    disp(sprintf('Monthly Average: file %d of %d complete',i,size(Files,1)));
end

amf = amf./cnt;
mno2 = mno2./cnt;
sno2 = sno2./cnt;
wcld = wcld./cnt;
err = err./cnt;
sc = sc./cnt;
if size(ave,3) > 6
    stratdata = stratdata./cnt;
end

v = version;
if str2num(v(1)) > 6
    if strcmp(Instrument, 'SCIAMACHY')
        weekmno2 = weekmno2./weekcnt;
        weeksno2 = weeksno2./weekcnt;
        %save([OutDir,ID,Instrument,num2str(gstep(1)),MonthName(Mn,:),num2str(Year),DirExt,'vert.mat'],'mno2','sno2','weekmno2','weeksno2','amf','wcld','err','cnt','sc','stratdata','lat','lon','-v6')
        save([OutDir,ID,Instrument,num2str(gstep(1)),MonthName(Mn,:),num2str(Year),DirExt,'vert',ExNeg,'.mat'],'mno2','sno2','weekmno2','weeksno2','amf','wcld','err','cnt','sc','lat','lon','-v6')
    elseif strcmp(Instrument, 'GOME')
        save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'vert',ExNeg,'.mat'],'mno2','sno2','amf','wcld','err','cnt','sc','stratdata','lat','lon','-v6')
    end
else
    if strcmp(Instrument, 'SCIAMACHY')
        weekmno2 = weekmno2./weekcnt;
        weeksno2 = weeksno2./weekcnt;
        %save([OutDir,ID,Instrument,num2str(gstep(1)),MonthName(Mn,:),num2str(Year),DirExt,'vert.mat'],'mno2','sno2','weekmno2','weeksno2','amf','wcld','err','sc','stratdata','lat','lon')
        save([OutDir,ID,Instrument,num2str(gstep(1)),MonthName(Mn,:),num2str(Year),DirExt,'vert',ExNeg,'.mat'],'mno2','sno2','weekmno2','weeksno2','amf','wcld','err','sc','lat','lon')
    elseif strcmp(Instrument, 'GOME')
        save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'vert',ExNeg,'.mat'],'mno2','sno2','amf','wcld','err','sc','stratdata','lat','lon')
    end    
end

if OutputHDF == 1
    hdfname = [OutDir,ID,Instrument,num2str(gstep(1)),MonthName(Mn,:),num2str(Year),DirExt,'vert',ExNeg,'.he5'];
    hdf5write(hdfname,'/lat',single(lat));        
    hdf5write(hdfname,'/lon',single(lon),'WriteMode','append');        
    hdf5write(hdfname,'/sno2',single(sno2),'WriteMode','append');        
    hdf5write(hdfname,'/mno2',single(mno2),'WriteMode','append');        
    hdf5write(hdfname,'/sc',single(sc),'WriteMode','append');        
end
