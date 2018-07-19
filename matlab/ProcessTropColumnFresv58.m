Species = 'NO2'; % NO2
Instrument = 'SCIAMACHY'; % SCIAMACHY / GOME
Source = 'Dal'; %Dal/KNMI

% Parent Input/Output Directories

% For Bastien's 2000 NO2
%InDir = '/data/bsauvage/AMF_fortran/amfv5.2.std/';
%OutDir = '/data/kelaar/amf/bsauvage/';

% Main GOME run
%InDir = '/raid-spang/nmoore/amf/supreme_amf_eta/';
%OutDir = '/raid-spang/kelaar/amf/GOME-20070716/';

% SCIAMACHY directories
%InDir = '/data1/kelaar/amf/amfv5.2.stetson/scia.v2.1.amf.simos/';
InDir = '/data3/akhila/in_progress_amfv5.8/SCIA-NO2/';
InDirSource = '/data1/kelaar/data/SCIAMACHY/Version2.1/';
%OutDir = '/data1/kelaar/amf/amfv5.2.stetson/scia.v2.1.simos/';
OutDir = '/data3/akhila/in_progress_amfv5.8/akhila_out/';

DirExt = ''; %extension onto date of subdirectories (input/output)
ID = '';

DailyStratDir = '/data1/kelaar/TEMIS/';
UseOMIstrat = 0;
wDiurnal = 0; % UseOMIstrat with diurnal variation
Windows = 0; % is running in windows environment? (changes directory '/' -> '\' for subdirectories)
Frescov5 = 1;
JustRegrid = 0; % only run the monthly regrid (load the rest)
newsmooth = 1;
OutputHDF = 0; % include monthly mean in hdf format
ExcludeNegatives = 0; % exclude negative values from monthly averages
amfv58 = 1;  % use amfv5.8 formatting (two separate input files)

MonthName = ['jan';'feb';'mar';'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'];

if strcmp(Species,'NO2')
    %FileExt = '.no2';
    FileExt = '.v5';
elseif strcmp(Species,'HCHO')
    FileExt = '.hcho.gomecat';
end

if ~strcmp(Instrument,'SCIAMACHY') && (amfv58 == 1) && ~strcmp(Source,'Dal')
    disp('amfv5.8 formatting only enabled for Dalhousie SCIAMACHY!')
    pause
end

for Year = [2005]
    for Mn = 3:3 % Month to process
		
        CarryOn = 1;
        
        if 0 == 1 %(Year == 2000 & Mn == 3) | (Year == 2000 & Mn == 4) | (Frescov5 == 1 & Year == 2005 & Mn == 10 | Year == 2005 & Mn == 3)
            FindValidFiles = 1; % inspect each file to ensure is valid (safety check for percieved change in # of columns, generally leave off)
        else
            FindValidFiles = 0;
        end
        %FindValidFiles = 0;
		
		if JustRegrid == 0
            if UseOMIstrat == 1
                if wDiurnal == 1
                    Yr = 2005
                    lokvstrat = hdfread(sprintf('/data1/lamsal/For_Aaron/OMI/bin_stdstrNO2_%d_new.hdf',Yr*100+Mn),'STRAT_NO2');
                    loklat = hdfread(sprintf('/data1/lamsal/For_Aaron/OMI/bin_stdstrNO2_%d_new.hdf',Yr*100+Mn),'LATITUDE');
                    loklon = hdfread(sprintf('/data1/lamsal/For_Aaron/OMI/bin_stdstrNO2_%d_new.hdf',Yr*100+Mn),'LONGITUDE');                    
                else
                    lokvstrat = hdfread(sprintf('/data1/lamsal/For_Aaron/bin_stdstrNO2_%d.hdf',Yr*100+Mn),'STRAT_NO2');
                    loklat = hdfread(sprintf('/data1/lamsal/For_Aaron/bin_stdstrNO2_%d.hdf',Yr*100+Mn),'LATITUDE');
                    loklon = hdfread(sprintf('/data1/lamsal/For_Aaron/bin_stdstrNO2_%d.hdf',Yr*100+Mn),'LONGITUDE');
                end
            end

            % Check if Slant Columns have already been calculated
            get_tslant_v2_new_Fresv58
            
        end

        if CarryOn == 1

            if JustRegrid == 0
                % Stage 2
                get_orbit_v2_new_Fresv58
            end

            if strcmp(Instrument, 'GOME')
                gstep = [2 2.5];
                %gstep = [0.5 0.5];
                get_month_v2
            elseif strcmp(Instrument, 'SCIAMACHY')
                gstep = [2 2.5];
                get_month_v2_Fres
                
                %gstep = [1 1.25];
                %get_month_v2_Fres

                gstep = [0.4 0.4];
                get_month_v2_Fres;
                
            end
            
        end
	end
end

disp('fin.')
