% Get tropospheric vertical column for each orbit

%mean orbital slant column
mnOrbit = median(Orbit(isfinite(Orbit)));

%stdev orbital slant column
stdOrbit = std(Orbit(isfinite(Orbit)));

gstep = 1;
lat = [-50:gstep:66];
lon = 0:gstep:360-gstep;

clear v
% format of output files for tropospheric vertical column
% SCIAMACHY columns are shifted relative to GOME by 'v'
% retrieve list of files
if strcmp(Instrument,'GOME')
    writeformat = '%5i %10.3e %10.3e %8.4f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.4f %10.3e %10.3e %10.3e\n';
    v=0;
    if Windows == 1
        Files = dir(sprintf('%s%d-%s%s\\*%s',InDir,Year,MonthName(Mn,:),DirExt,FileExt));
    else
        Files = dir(sprintf('%s%d-%s%s/*%s',InDir,Year,MonthName(Mn,:),DirExt,FileExt));
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
    
    if strcmp(Species,'NO2')
        MaxDSC = 4.5e15;
    else
        MaxDSC = 1e18;
    end
elseif strcmp(Instrument,'SCIAMACHY')
    if Frescov5 == 0
        writeformat = '%5i %10.3e %10.3e %8.4f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.4f %10.3e %10.3e %10.3f %5i %10.3e\n';
    else
        writeformat = '%5i %10.3e %10.3e %8.4f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.4f %10.3e %10.3e %10.3f %5i %10.3e\n';        
    end
    v=1;
    if strcmp(Source,'Dal')
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
    elseif strcmp(Source,'KNMI')
        if Windows == 1
            Files = dir(sprintf('%s%d-%s%s\\SCIA-*',InDir,Year,MonthName(Mn,:),DirExt));
        else
            Files = dir(sprintf('%s%d-%s%s/SCIA-*',InDir,Year,MonthName(Mn,:),DirExt));
        end        
    end
    MaxDSC = 3e15;
end

% Create Output Directory (if needed)
if Windows == 1
    temp = dir(sprintf('%s%d-%s%s\\',OutDir,Year,MonthName(Mn,:),DirExt));
else
    temp = dir(sprintf('%s%d-%s%s/',OutDir,Year,MonthName(Mn,:),DirExt));
end
if isempty(temp)
    mkdir(OutDir,sprintf('%d-%s%s',Year,MonthName(Mn,:),DirExt));
end

for i = 1:size(Files,1)
    if ~FindValidFiles
        if amfv58 ~= 1
            if Windows == 1
                dat = load(sprintf('%s%d-%s%s\\%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
            else
                dat = load(sprintf('%s%d-%s%s/%s',InDir,Year,MonthName(Mn,:),DirExt,Files(i).name));
            end
        else
            % Combine amf output with original files to mirror old
            % formatting
            [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16] = textread(sprintf('%s/%d-%s%s/%s',InDir,Year,strMn(end-1:end),DirExt,Files(i).name),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');
            [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24] = textread(sprintf('%s%d-%s%s/%s',InDirSource,Year,MonthName(Mn,:),DirExt,Files(i).name(5:end-4)),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]','delimiter',' ');

            A = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16];
            B = [B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24];
            [C, IA, IB] = intersect(A1*1000+A2,B1*1000+B2);
            dat = [A(IA,1:2) zeros(size(IA)) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(IA)) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
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
                dat = str2num(fline);
                numcol = length(Data);
                while(~feof(fid))
                    fline = fgetl(fid);
                    temp = str2num(fline);
                    if numcol == length(temp)
                        dat = [Data; temp];
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
            [C, IA, IB] = intersect(A(1,:)*1000+A(2,:),B(1,:)*1000+B(2,:));
            dat = [A(IA,1:2) zeros(size(A,1),1) B(IB,21:23) A(IA,3:4) B(IB,16:20) zeros(size(A,1),1) A(IA,7:8) B(IB,6:15) A(IA,11:13) A(IA,15:16)];
            clear A B temp temp2
        end
    end
            
    if (length(dat)>0)
        if Frescov5 == 0
            dat = dat((dat(:,21+v)<65 & dat(:,21+v)>-49),:); %filter for latitude
        else
            dat = dat((dat(:,24+v)<65 & dat(:,24+v)>-49),:); %filter for latitude            
        end
    end

    if strcmp(Instrument,'SCIAMACHY')
        proceed = (length(dat)>0 & OrbitStat.Filter(i)~=1);
    elseif strcmp(Instrument,'GOME')
        GOME_FileSpot = GOME_FileSpot+1;
        if strcmp(Files(i).name,GOME_Filter(GOME_FileSpot).FileName)
            proceed = (length(dat)>0 & (GOME_Filter(GOME_FileSpot).Valid==1));
        else
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
            else
                proceed = (length(dat)>0 & (GOME_Filter(GOME_FileSpot).Valid==1));
            end
        end
    end
    
    if proceed
        if Windows == 1
            fname = sprintf('%s%d-%s%s\\%s',OutDir,Year,MonthName(Mn,:),DirExt,Files(i).name);
        else
            fname = sprintf('%s%d-%s%s/%s',OutDir,Year,MonthName(Mn,:),DirExt,Files(i).name);
        end
        fid = fopen(fname,'w');
        for j = 1:size(dat,1)

            if Frescov5 == 0
                [temp latind] = min(abs(lat-dat(j,21+v)));
                [temp lonind] = min(abs(lon-dat(j,22+v)));
            else
                [temp latind] = min(abs(lat-dat(j,24+v)));
                [temp lonind] = min(abs(lon-dat(j,25+v)));                
            end
			gamf = sec(pi/180*dat(j,6+v)) + sec(pi/180*dat(j,7+v));
            if strcmp(Instrument,'SCIAMACHY')
                if OrbitStat.Shifted(i) == 1
                    dat(j,3+v) = dat(j,3+v) - OrbitStat.Shift;
                end
            elseif strcmp(Instrument,'GOME')
                
                dat(j,3+v) = dat(j,3+v) + GOME_Filter(GOME_FileSpot).Shift;
                
            end
            
            SC = dat(j,3+v);
            if UseOMIstrat == 0
                if Frescov5 == 0
                    trop = (dat(j,3+v) - strat(latind)*strat_scale_2d(latind,lonind)*gamf) / dat(j,25+v); % Tropospheric NO2 Column [molec/cm2]
                    stratslant = strat(latind)*strat_scale_2d(latind,lonind)*gamf;
                else
                    trop = (dat(j,3+v) - strat(latind)*strat_scale_2d(latind,lonind)*gamf) / dat(j,28+v); % Tropospheric NO2 Column [molec/cm2] 
                    stratslant = strat(latind)*strat_scale_2d(latind,lonind)*gamf;
                end
            elseif UseOMIstrat == 1
                if Frescov5 == 0
                    trop = (dat(j,3+v) - OMIstrat(latind,lonind)*gamf) / dat(j,25+v); % Tropospheric NO2 Column [molec/cm2]
                    stratslant = OMIstrat(latind,lonind)*gamf;
                else
                    trop = (dat(j,3+v) - OMIstrat(latind,lonind)*gamf) / dat(j,28+v); % Tropospheric NO2 Column [molec/cm2]  
                    stratslant = OMIstrat(latind,lonind)*gamf;
                end                
            end
            err = sqrt((dat(j,4+v)/gamf)^2 + (3e14)^2 + (0.4*trop)^2); % Total uncertainty [molec/cm2] 
    
            Pixel = dat(j,1);
            scan = dat(j,2);
            DSC = dat(j,4+v); %deviation of slant column
            SC = dat(j,3+v); %slant column
            if Frescov5 == 0
                FractionalRadianceFromClouds = dat(j,26+v);
                AMF = dat(j,25+v);
                GEOS_VC = dat(j,27+v);

                Lat1 = dat(j,13+v);
                Lon1 = dat(j,14+v);
                Lat2 = dat(j,15+v);
                Lon2 = dat(j,16+v);
                Lat3 = dat(j,17+v);
                Lon3 = dat(j,18+v);
                Lat4 = dat(j,19+v);
                Lon4 = dat(j,20+v);
                LatCentre = dat(j,21+v);
                LonCentre = dat(j,22+v);
                
            else
                FractionalRadianceFromClouds = dat(j,29+v);
                AMF = dat(j,28+v);
                GEOS_VC = dat(j,30+v);

                Lat1 = dat(j,16+v);
                Lon1 = dat(j,17+v);
                Lat2 = dat(j,18+v);
                Lon2 = dat(j,19+v);
                Lat3 = dat(j,20+v);
                Lon3 = dat(j,21+v);
                Lat4 = dat(j,22+v);
                Lon4 = dat(j,23+v);
                LatCentre = dat(j,24+v);
                LonCentre = dat(j,25+v);
                
            end
            
            % remove backscan and bad pixels
            if (Lon1 > Lon4) & (DSC < MaxDSC)
                if strcmp(Instrument,'SCIAMACHY')
                    UTC = dat(j,3); %UTC Time HHMMSS.### where ### is fraction of second
					fprintf(fid,writeformat,Pixel,SC,err,FractionalRadianceFromClouds,Lat1,Lon1,Lat2,Lon2,Lat3,Lon3,Lat4,Lon4,LatCentre,LonCentre,AMF,GEOS_VC,trop,UTC,scan,stratslant);
                elseif strcmp(Instrument, 'GOME')
					fprintf(fid,writeformat,Pixel,SC,err,FractionalRadianceFromClouds,Lat1,Lon1,Lat2,Lon2,Lat3,Lon3,Lat4,Lon4,LatCentre,LonCentre,AMF,GEOS_VC,trop,stratslant);
                end
            end
        end
        fclose(fid);
    end
    disp(sprintf('Get Orbit: file %d of %d complete',i,size(Files,1)));    
end
