if strcmp(Instrument,'GOME')
    if Windows == 1
        Files = dir(sprintf('%s%d-%s%s\\*%s',InDir,Year,MonthName(Mn,:),DirExt,FileExt));
    else
        Files = dir(sprintf('%s%d-%s%s/*%s',InDir,Year,MonthName(Mn,:),DirExt,FileExt));
    end
elseif strcmp(Instrument,'SCIAMACHY')
    if strcmp(Source,'KNMI')
        if Windows == 1
            Files = dir(sprintf('%s%d-%s%s\\SCIA-*',InDir,Year,MonthName(Mn,:),DirExt));
        else
            Files = dir(sprintf('%s%d-%s%s/SCIA-*',InDir,Year,MonthName(Mn,:),DirExt));
        end        
    elseif strcmp(Source,'Dal')
        if amfv58 == 1
            strMn = num2str(Mn+1000);
            Files = dir(sprintf('%s/%d-%s%s/AMF-orb*%s.no2',InDir,Year,strMn(end-1:end),DirExt,FileExt));
        else
            if Windows == 1
                Files = dir(sprintf('%s%d-%s%s\\orb*',InDir,Year,MonthName(Mn,:),DirExt));
            else
                Files = dir(sprintf('%s%d-%s%s/orb*',InDir,Year,MonthName(Mn,:),DirExt));
            end
        end
    end
end

if isempty(Files)
    CarryOn = 0; 
end

if CarryOn == 1

    gstep=1;
    Lat = [-50:gstep:66];
    Lon = 0:gstep:361-gstep;

    cnt = zeros(length(Lat),length(Lon));
    no2 = zeros(length(Lat),length(Lon));
    
    if UseOMIstrat == 1
        latspot = [find(Lat(1) == loklat) find(Lat(end) == loklat)] ;
        OMIstrat = [lokvstrat(latspot(1):latspot(end),:) lokvstrat(latspot(1):latspot(end),1)];
    elseif UseOMIstrat == 0
        vstrat = zeros(length(Lat),length(Lon));        
    end

    % SCIAMACHY columns are shifted relative to GOME by 'v'
    clear v
    if strcmp(Instrument,'SCIAMACHY')
        v = 1;
    elseif strcmp(Instrument,'GOME')
        v = 0;
    end

    if strcmp(Instrument,'SCIAMACHY')
        for i = 1:size(Files,1)

            if ~FindValidFiles
                if amfv58 ~= 1
                    if Windows == 1
                        Data = load(sprintf('%s%d-%s%s\\%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                    else
                        Data = load(sprintf('%s%d-%s%s/%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                    end
                else
                    % Combine amf output with original files to mirror old
                    % formatting
                    test=(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name))
                    [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16] = textread(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');
                    test2=(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)))
                    [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24] = textread(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');
                    
                    A = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16];
                    B = [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24];
                    [C, IA, IB] = intersect(A1*1000+A2,B1*1000+B2);
                    Data = [A(IA,1:2) zeros(size(IA)) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(IA)) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
                    clear A A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16
                    clear B B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24
                end
            else
                if amfv58 ~= 1
                    if Windows == 1
                        fid = fopen(sprintf('%s%d-%s%s\\%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));            
                    else
                        fid = fopen(sprintf('%s%d-%s%s/%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                    end
                    fline = fgetl(fid);
                    if fline ~= -1
                        Data = str2num(fline);
                        numcol = length(Data);
                        while(~feof(fid))
                            fline = fgetl(fid);
                            temp = str2num(fline);
                            if numcol == length(temp)
                                Data = [Data; temp];
                            end
                        end
                    end
                    fclose(fid);        
                else
                    for v58i = 1:2
                        if v58i == 1
                            fid = fopen(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name));
                        else
                            fid = fopen(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)));
                        end
                        fline = fgetl(fid);
                        
                        if fline ~= -1
                        
                            temp = sscanf(fline,'%f','delimiter',' ');

                            numcol = length(temp);
                            while(~feof(fid))
                                fline = fgetl(fid);
                                fline2 = fgetl(fid2);
                                temp2 = sscanf(fline,'%f','delimiter',' ');
                                if numcol == length(A)
                                    temp = [temp; temp2];
                                end
                            end
                        end
                        fclose(fid);    
                        
                        if v58i == 1
                            A = temp;
                        else
                            B = temp;
                        end
                    end
                    [C, IA, IB] = intersect(A(1,:)*1000+A(2,:),B(2,:)*1000+B(2,:));
                    Data = [A(IA,1:2) zeros(size(A,1),1) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(A,1),1) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
                    clear A B temp temp2
                end
            end

            if ~isempty(Data)
                % Filter points outside of specified latitude range
                if Frescov5 == 0
                    Data = Data((Data(:,21+v) < (max(Lat) - gstep) & Data(:,21+v) > (min(Lat) + gstep)),:);
                else
                    Data = Data((Data(:,24+v) < (max(Lat) - gstep) & Data(:,21+v) > (min(Lat) + gstep)),:);                    
                end
            end

            if ~isempty(Data)

                Orbit(i) = mean(Data(:,3+v));
                Error(i) = mean(Data(:,4+v));

                if Error(i) > 4.5e15 & strcmp(Instrument,'SCIAMACHY')
                    Error(i) = NaN;
                end
            else
                Orbit(i) = NaN;
                Error(i) = NaN;
            end

            disp(sprintf('Get Filter: file %d of %d complete',i,size(Files,1)));

        end

        orb = Orbit;

        % find plateaus
        k = 1;
        j = 10;
        c = 0;
        spot = [1];
        while j<length(orb)
            ave = nanmean(orb(k:j));
            stdev = nanstd(orb(k:j));
            if ~isnan(orb(j))
                if ave+stdev > orb(j) & ave-stdev < orb(j)
                    c = 0;
                else
                    c = c + 1;
                end
            end

            j = j+1;

            if c == 10
                spot = [spot; j-10];
                k = j-10;
                j = k+10;
                c = 0;
            end
        end

        forb = orb;
        OrbitStat.Start = spot;
        OrbitStat.Shifted = zeros(size(orb));
        OrbitStat.Ave = [];
        OrbitStat.Stdev = [];

        if length(spot)>1
            % determine which state shows heightened levels
            tempspot1 = [];
            for j = 1:2:length(spot)
                if j < length(spot)
                    tempspot1 = [tempspot1 [spot(j):spot(j+1)-1]];
                else
                    tempspot1 = [tempspot1 [spot(j):length(orb)]];
                end
            end

            tempspot2 = [];
            for j = 2:2:length(spot)
                if j < length(spot)
                    tempspot2 = [tempspot2 [spot(j):spot(j+1)-1]];
                else
                    tempspot2 = [tempspot2 [spot(j):length(orb)]];
                end
            end

            OrbitStat.Shift = abs((nanmedian(orb(tempspot2)) - nanmedian(orb(tempspot1))));
            if nanmedian(orb(tempspot2)) > nanmedian(orb(tempspot1))
                On = 2;
                forb(tempspot2) = forb(tempspot2) - OrbitStat.Shift;
                OrbitStat.Shifted(tempspot2) = 1;
            else
                On = 1;
                forb(tempspot1) = forb(tempspot1) - OrbitStat.Shift;
                OrbitStat.Shifted(tempspot1) = 1;
            end

        end

        OrbTemp.data = [];
        temp = forb;
        oldstd = -1;
        while(1)
            tempave = nanmean(temp);
            tempstd = nanstd(temp);
            temp(temp>tempave+3*tempstd|temp<tempave-3*tempstd) = NaN;
            OrbTemp.data = [OrbTemp.data; [tempave tempstd]];
            if oldstd == tempstd
                break
            end
            oldstd = tempstd;
        end

        OrbitStat.Ave = [OrbitStat.Ave; tempave];
        OrbitStat.Stdev = [OrbitStat.Stdev; tempstd];
        OrbitStat.Filter = isnan(temp);

        forb = temp;

    elseif strcmp(Instrument,'GOME')

        if strcmp(Source,'Dal')
            load(sprintf('GOME-Filter-%s',Species));
        elseif strcmp(Source,'KNMI')
            load(sprintf('%s-GOME-Filter-%s',Source,Species));        
        end

        % Find File Location in GOME_Filter
        GOME_FileSpot = NaN;
        for i = 1:length(GOME_Filter)
            if strcmp(GOME_Filter(i).FileName,Files(1).name)
                GOME_FileSpot = i;
                break
            end
        end

        if isnan(GOME_FileSpot)
            disp('GOME Filter Error: File not found in filter list!');
            pause
        end    

    end

    for i = 1:size(Files,1)
        if ~FindValidFiles
            if amfv58 ~= 1
                if Windows == 1
                    Data = load(sprintf('%s%d-%s%s\\%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                else
                    Data = load(sprintf('%s%d-%s%s/%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                end
            else
                % Combine amf output with original files to mirror old
                % formatting
                [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16] = textread(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');
                [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24] = textread(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');

                A = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16];
                B = [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24];
                [C, IA, IB] = intersect(A1*1000+A2,B1*1000+B2);
                Data = [A(IA,1:2) zeros(size(IA)) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(IA)) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
                %disp(sprintf('i = %d, size(Data) = %d, %d',i,size(Data,1),size(Data,2)));
                clear A A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16
                clear B B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24
            end
        else
            if amfv58 ~= 1
                if Windows == 1
                    fid = fopen(sprintf('%s%d-%s%s\\%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));            
                else
                    fid = fopen(sprintf('%s%d-%s%s/%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
                end
                fline = fgetl(fid);
                if fline ~= -1
                    Data = str2num(fline);
                    numcol = length(Data);
                    while(~feof(fid))
                        fline = fgetl(fid);
                        temp = str2num(fline);
                        if numcol == length(temp)
                            Data = [Data; temp];
                        end
                    end
                end
                fclose(fid);        
            else
                for v58i = 1:2
                    if v58i == 1
                        fid = fopen(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name));
                    else
                        fid = fopen(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)));
                    end
                    fline = fgetl(fid);

                    if fline ~= -1

                        temp = sscanf(fline,'%f','delimiter',' ');

                        numcol = length(temp);
                        while(~feof(fid))
                            fline = fgetl(fid);
                            fline2 = fgetl(fid2);
                            temp2 = sscanf(fline,'%f','delimiter',' ');
                            if numcol == length(A)
                                temp = [temp; temp2];
                            end
                        end
                    end
                    fclose(fid);    

                    if v58i == 1
                        A = temp;
                    else
                        B = temp;
                    end
                end
                [C, IA, IB] = intersect(A(1,:)*1000+A(2,:),B(2,:)*1000+B(2,:));
                Data = [A(IA,1:2) zeros(size(IA)) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(IA)) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
                clear A B temp temp2
            end
        end

        if ~isempty(Data)
            % Filter points outside of specified latitude range
            if Frescov5 == 0
                Data = Data((Data(:,21+v) < (max(Lat) - gstep) & Data(:,21+v) > (min(Lat) + gstep)),:);
            else
                Data = Data((Data(:,24+v) < (max(Lat) - gstep) & Data(:,24+v) > (min(Lat) + gstep)),:);                
            end
        end

        % check GOME file in filter list
        if strcmp(Instrument,'GOME')
            GOME_FileSpot = GOME_FileSpot+1;
            if ~strcmp(Files(i).name,GOME_Filter(GOME_FileSpot).FileName)
                GOME_FileSpot = NaN;
                for q = 1:length(GOME_Filter)
                    if strcmp(GOME_Filter(q).FileName,Files(i).name)
                        GOME_FileSpot = q;
                        break
                    end
                end

                if isnan(GOME_FileSpot)
                    disp('GOME Filter File Mismatch!')
                    pause
                end
            end
        end

        if ~isempty(Data)

            Orbit(i) = mean(Data(:,3+v));
            Error(i) = mean(Data(:,4+v));

            if Error(i) > 4.5e15 & strcmp(Instrument,'SCIAMACHY')
                Error(i) = NaN;
            end


            for j = 1:size(Data,1)

                if strcmp(Instrument,'SCIAMACHY')
                    if OrbitStat.Filter(i) ~= 1

                        if Frescov5 == 0
                            iLat = floor(Data(j,21+v)-min(Lat)+1);
                            iLon = floor(Data(j,22+v)-min(Lon)+1);
                        else
                            iLat = floor(Data(j,24+v)-min(Lat)+1);
                            iLon = floor(Data(j,25+v)-min(Lon)+1);                            
                        end
                        cnt(iLat,iLon) = cnt(iLat,iLon) + 1;
                        gamf = sec(pi/180*Data(j,6+v)) + sec(pi/180*Data(j,7+v));

                        if OrbitStat.Shifted(i) == 1
                            Data(j,3+v) = Data(j,3+v) - OrbitStat.Shift;
                        end

                        no2(iLat,iLon) = no2(iLat,iLon) + Data(j,3+v);  
                        if UseOMIstrat == 0
                            if Frescov5 == 0
                                vstrat(iLat,iLon) = vstrat(iLat,iLon) + (Data(j,3+v)-Data(j,25+v)*Data(j,27+v))/gamf;
                            else
                                % vstrat = vstrat + (SCIA Slant Column - Model Slant Column) / AMF
                                % vstrat = vertical stratospheric column
                                vstrat(iLat,iLon) = vstrat(iLat,iLon) + (Data(j,3+v)-Data(j,28+v)*Data(j,30+v))/gamf;                            
                            end
                        end
                    end
                elseif strcmp(Instrument,'GOME')

                    % Shift Data for GOME Filter
                    Data(j,3+v) = Data(j,3+v) + GOME_Filter(GOME_FileSpot).Shift;

                    iLat = floor(Data(j,21+v)-min(Lat)+1);
                    iLon = floor(Data(j,22+v)-min(Lon)+1);
                    cnt(iLat,iLon) = cnt(iLat,iLon) + 1;
                    no2(iLat,iLon) = no2(iLat,iLon) + Data(j,3+v);  
                    gamf = sec(pi/180*Data(j,6+v)) + sec(pi/180*Data(j,7+v));
                    
                    if UseOMIstrat == 0
                        vstrat(iLat,iLon) = vstrat(iLat,iLon) + (Data(j,3+v)-Data(j,25+v)*Data(j,27+v))/gamf;
                    end
                end
            end

        else
            Orbit(i) = NaN;
            Error(i) = NaN;
        end

        disp(sprintf('Get TSlant: file %d of %d complete',i,size(Files,1)));

    end

    % Average values
    cnt(cnt==0)=nan;
    tslant = no2./cnt;
    
    if UseOMIstrat == 0
        vstrat = vstrat./cnt;

        if newsmooth == 1
            strat = get_strat_1d_new(Lat,Lon,vstrat);
        else
            strat = get_strat_1d(Lat,Lon,vstrat);
        end
        get_strat_2d;
        OMIstrat = NaN;
    elseif UseOMIstrat == 1
        strat_scale_2d = NaN;
        vstrat = NaN;
        strat = NaN;
    end

    v = version;
    if str2num(v(1)) > 6
        if strcmp(Instrument,'SCIAMACHY')
            save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'slant'],'vstrat','tslant','strat','Orbit','Error','OrbitStat','strat_scale_2d','OMIstrat','-v6')
        elseif strcmp(Instrument,'GOME')
            save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'slant'],'vstrat','tslant','strat','Orbit','Error','GOME_Filter','OMIstrat','-v6')
        end
    else
        if strcmp(Instrument,'SCIAMACHY')
            save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'slant'],'vstrat','tslant','strat','Orbit','Error','OrbitStat','strat_scale_2d')
        elseif strcmp(Instrument,'GOME')
            save([OutDir,ID,Instrument,MonthName(Mn,:),num2str(Year),DirExt,'slant'],'vstrat','tslant','strat','Orbit','Error','GOME_Filter')
        end
    end
    disp('fin.')
end