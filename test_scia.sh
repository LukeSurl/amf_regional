#!/bin/csh -f

#
# Script to loop over all the GOME slant column files.
#
limit stacksize 976m



set month = {"jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"}

set days = {"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"}

#set potfiles = {"OMI-Aura_L2-OMHCHO_2012m0420t0104-o41296_v003-2014m0701t025437.txt","OMI-Aura_L2-OMHCHO_2012m0420t0243-o41297_v003-2014m0701t025920.txt","OMI-Aura_L2-OMHCHO_2012m0420t0421-o41298_v003-2014m0701t025211.txt","OMI-Aura_L2-OMHCHO_2012m0420t0600-o41299_v003-2014m0701t025100.txt","OMI-Aura_L2-OMHCHO_2012m0420t0739-o41300_v003-2014m0701t025646.txt","OMI-Aura_L2-OMHCHO_2012m0420t0918-o41301_v003-2014m0701t030703.txt","OMI-Aura_L2-OMHCHO_2012m0420t1057-o41302_v003-2014m0701t030214.txt","OMI-Aura_L2-OMHCHO_2012m0420t1236-o41303_v003-2014m0701t030532.txt","OMI-Aura_L2-OMHCHO_2012m0420t1415-o41304_v003-2014m0701t030223.txt","OMI-Aura_L2-OMHCHO_2012m0420t1554-o41305_v003-2014m0701t030436.txt","OMI-Aura_L2-OMHCHO_2012m0420t1733-o41306_v003-2014m0701t031104.txt","OMI-Aura_L2-OMHCHO_2012m0420t1911-o41307_v003-2014m0701t030205.txt","OMI-Aura_L2-OMHCHO_2012m0420t2050-o41308_v003-2014m0701t030255.txt","OMI-Aura_L2-OMHCHO_2012m0420t2229-o41309_v003-2014m0701t030615.txt"}

set DIR = "/group_workspaces/jasmin/geoschem/local_users/lsurl/OMI_ASCII/CH2O"
set Output = "/home/users/lsurl/amf581p/output"

set yr =  "20"$days[$argv[1]]      #year
set mon = $days[$argv[2]]    #two digit month (i.e. 07 for July)
set daycount = $argv[3]      #first day to be run
set endDay = $argv[4]   #last day to be run
#set thisfile = $potfiles[$argv[5]]

set molecule = 0     #0=HCHO, 1=NO2
set append = '.bpch'
set FRESCO5 = 0   #1= FRESCO v5, 0= old Fresco

#set prefix = '/data3/akhila/Code/GEOS-Chem.v8-02-02.stdrun/run.v8-02-02/SCIA_runs/ts_satellite_9_12.'
set prefix = '/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_2x25_fullchem/ts_satelliteOMI.'

#set tropfile = '/data3/akhila/Code/GEOS-Chem.v8-02-02.stdrun/run.v8-02-02/trophghts.20050101.bpch'
set tropfile = '/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_2x25_fullchem/tp-hght-201301'
#set tropfile = '/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/TPheight_05.bpch'
#set DIR = $DIR$yr-$month[$mon] #Don't change dir here 

while ( $daycount <= $endDay )

set infile = $prefix$yr$mon$days[$daycount]$append
set aerfile = $prefix$yr$mon$days[$daycount]$append
set cloudFile = 'dummyFile'

#echo "FILE:"
#echo $DIR/orbSTAR-$yr$mon$days[$daycount]

#echo $DIR/$yr/$yr$mon/"GOME2_MetOpA_UoL_HCHO_"$yr$mon$days[$daycount]
#foreach file ( $DIR/$yr/$yr$mon/$thisfile ) #Just going to have to scan every file in folder-no way to cleanly pick out days using filenames 
foreach file ($DIR/$yr/$yr$mon/* )

#skip if file not referring to this day

set fileday = `echo $file | awk '{print substr ($0,104,4) }'` 

if ( "$fileday" != "$mon$days[$daycount]") then
 echo $fileday
 continue
endif 


#skip if file is empty
if (! -s $file ) then
 continue
endif

#skip if an aux file
if ( "`echo $file | cut -d'.' -f2`" == "aux" )  then
  echo $file is a aux file
  continue
else
  echo $file is not a zip file
endif




set filetag = ` echo $file | sed 's,'$DIR/',,'`
#set filetag = $Output"/"$yr"-"$mon"/AMF-"$filetag".HCHO"
set filetag = $Output"/AMF/"$filetag".hcho"
#set outfile = $Output"/"$yr"-"$mon"/"
set outfile = $Output"/"
set outputFile_wscat = ` echo $filetag | sed 's,'.hcho','.wscat','`
set outputFile_shape = ` echo $filetag | sed 's,'.hcho','.shape','`

#echo "FILETAG:"
echo $filetag
#echo "OUTPUTFILE_WSCAT;"
echo $outputFile_wscat
#echo "OUTPUTFILE_SHAPE"
echo $outputFile_shape
#echo "INFILE"
echo $infile
#echo "TROPFILE"
echo $tropfile

./amf.run<<EOF
$file
$cloudFile
$days[$daycount]
$infile
$aerfile
$tropfile
$filetag
$outputFile_wscat
$outputFile_shape
$yr
$mon
$molecule
$FRESCO5
EOF


end # for 

@ daycount++
end # while

