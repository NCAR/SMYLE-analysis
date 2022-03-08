# trak TCs from SMYLE hindcasts, example for Feb starts

#!/bin/bash -l

###=======================================================================
#PBS -N tempest.par
#PBS -A P93300313 
#PBS -l walltime=00:05:00
#PBS -q regular
#PBS -j oe
#PBS -l select=1:ncpus=36:mpiprocs=36
###############################################################

############ USER OPTIONS #####################

## Unique string (useful for processing multiple data sets in same folder
SIMYRS=${TRACK_YEAR}     # 1979_2012, RCP85_2070_2099,
SIMYRSEND=`expr $SIMYRS + 2`
ENSMEM=${ENS_MEM}
INIMONTH=02

echo ${SIMYRS}
echo ${SIMYRSEND}
echo ${ENSMEM}
## Path to TempestExtremes binaries on YS
TEMPESTEXTREMESDIR=/glade/work/zarzycki/tempestextremes/

## Topography filter file (needs to be on same grid as PSL, U, V, etc. data
TOPOFILE=/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_sgh30_24km_GRNL_c170103.nc

## If using unstructured CAM-SE ne120 data
CONNECTFLAG="" 


## Path where files are
PATHTOFILES=/glade/scratch/huili7/SMYLE_data


############ TRACKER MECHANICS #####################
starttime=$(date -u +"%s")

DATESTRING=`date +"%s%N"`
FILELISTNAME=filelist.txt.${DATESTRING}
touch ./$FILELISTNAME
TRAJFILENAME=trajectories.txt.SMYLE.${SIMYRS}.${INIMONTH}.0${ENSMEM}


# for feb start
# for feb start
PSLFILE=${PATHTOFILES}/b.e21.BSMYLE.f09_g17.${SIMYRS}-${INIMONTH}.0${ENSMEM}.cam.h2.PSL.${SIMYRS}020100-${SIMYRSEND}013100.nc
UBOTFILE=${PATHTOFILES}/b.e21.BSMYLE.f09_g17.${SIMYRS}-${INIMONTH}.0${ENSMEM}.cam.h2.UBOT.${SIMYRS}020100-${SIMYRSEND}013100.nc
VBOTFILE=${PATHTOFILES}/b.e21.BSMYLE.f09_g17.${SIMYRS}-${INIMONTH}.0${ENSMEM}.cam.h2.VBOT.${SIMYRS}020100-${SIMYRSEND}013100.nc
Z500FILE=${PATHTOFILES}/b.e21.BSMYLE.f09_g17.${SIMYRS}-${INIMONTH}.0${ENSMEM}.cam.h2.Z500.${SIMYRS}020100-${SIMYRSEND}013100.nc
Z300FILE=${PATHTOFILES}/b.e21.BSMYLE.f09_g17.${SIMYRS}-${INIMONTH}.0${ENSMEM}.cam.h2.Z300.${SIMYRS}020100-${SIMYRSEND}013100.nc
echo "${PSLFILE};${UBOTFILE};${VBOTFILE};${Z500FILE};${Z300FILE};${TOPOFILE}" >> $FILELISTNAME


DCU_PSLFOMAG=200.0
DCU_PSLFODIST=8
DCU_WCFOMAG=-6.    # Z300Z500 -6.0, T400 -0.4
DCU_WCFODIST=8.
DCU_WCMAXOFFSET=3.0
DCU_WCVAR="_DIFF(Z300,Z500)"   #DCU_WCVAR generally _DIFF(Z300,Z500) or T400
DCU_MERGEDIST=6.0
SN_TRAJRANGE=10.0
SN_TRAJMINLENGTH=10
SN_TRAJMAXGAP=3
SN_MAXTOPO=150.0
SN_MAXLAT=50.0
SN_MINWIND=8.0
SN_MINLEN=9



STRDETECT="--verbosity 0 --timestride 1 ${CONNECTFLAG} --out cyclones_tempest.${DATESTRING} --closedcontourcmd PSL,${DCU_PSLFOMAG},${DCU_PSLFODIST},0;${DCU_WCVAR},${DCU_WCFOMAG},${DCU_WCFODIST},${DCU_WCMAXOFFSET} --mergedist ${DCU_MERGEDIST} --searchbymin PSL --outputcmd PSL,min,0;_VECMAG(UBOT,VBOT),max,2;PHIS,max,0"
echo $STRDETECT

#touch cyclones.${DATESTRING}
mpiexec_mpt ${TEMPESTEXTREMESDIR}/bin/DetectNodes --in_data_list "${FILELISTNAME}" ${STRDETECT} </dev/null
cat cyclones_tempest.${DATESTRING}* >> cyclones.${DATESTRING}
rm cyclones_tempest.${DATESTRING}*

# stitch node for each of the cyclone candidate file

# Stitch candidate cyclones together
${TEMPESTEXTREMESDIR}/bin/StitchNodes --format "i,j,lon,lat,slp,wind,phis" --range ${SN_TRAJRANGE} --minlength ${SN_TRAJMINLENGTH} --maxgap ${SN_TRAJMAXGAP} --in cyclones.${DATESTRING} --out ${TRAJFILENAME} --threshold "wind,>=,${SN_MINWIND},${SN_MINLEN};lat,<=,${SN_MAXLAT},${SN_MINLEN};lat,>=,-${SN_MAXLAT},${SN_MINLEN};phis,<=,${SN_MAXTOPO},${SN_MINLEN}"


rm ${FILELISTNAME}
rm log*.txt
rm cyclones.${DATESTRING}

endtime=$(date -u +"%s")
tottime=$(($endtime-$starttime))
printf "${tottime},${TRAJFILENAME}\n" >> timing.txt






