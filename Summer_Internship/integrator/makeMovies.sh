#!/bin/bash

# ------------------------------------------ #
#                                            #
#               makeMovies.sh                #
#                                            #
# ------------------------------------------ #

MATLAB_CMD="matlab"
INTEGRATOR_CMD="integrator.out"
PATH_TO_SCRIPT=/home/miles/MATLAB/Summer_2012/
PATH_TO_INTEGRATOR=/home/miles/Summer_2012/Summer_Internship/integrator/
PATH_TO_OUTPUT=${PATH_TO_INTEGRATOR}output/
PATH_TO_MOVIES=${PATH_TO_INTEGRATOR}movies/
CWD=`pwd`

P=
E=
S=
Z=

while [ $# -gt 0 ]
do
    case "$1" in
        -p)  P=$2; shift;;
	    -e)  E=$2; shift;;
	    -z)  Z=$2; shift;;
        -r)  R=$2; shift;;
        *) break;;
    esac
    shift
done

if [[ -z $P || -z $E || -z $R || -z $Z ]]
then
    echo -n "Need to specify probability (-p), strain (-e), rate (-r), and size"
    echo " (-z)."
    exit 1
else

    # --- Make the executable

    make -C ${PATH_TO_INTEGRATOR}
    make clean -C ${PATH_TO_INTEGRATOR}
    
    if [[ $? == "0" ]]
    then
        
        mkdir -p ${PATH_TO_MOVIES}

        # --- Run integrator script with specified parameters

        ${PATH_TO_INTEGRATOR}${INTEGRATOR_CMD} -p $P -e $E -r $R -z $Z

        cd ${PATH_TO_OUTPUT}

        # --- Generate images with matlab and make into movie

        ${MATLAB_CMD} < ${PATH_TO_SCRIPT}plotNodes.m
        ls | grep .tiff | sort -n > movie.txt
        mencoder mf://@movie.txt -ovc x264 -vf scale=1024:768 -mf type=tiff -o ${PATH_TO_MOVIES}network_movie_z${Z/./_}p${P/./_}e${E/./_}r${R/./_}.avi

        ${MATLAB_CMD} < ${PATH_TO_SCRIPT}plotNonAffine.m
        ls | grep .tiff | sort -n > movie.txt
        mencoder mf://@movie.txt -ovc x264 -vf scale=1024:768 -mf type=tiff -o ${PATH_TO_MOVIES}nonaff_motion_z${Z/./_}p${P/./_}e${E/./_}r${R/./_}.avi
        
        ${MATLAB_CMD} < ${PATH_TO_SCRIPT}affine_and_nonaff.m
        ls | grep .tiff | sort -n > movie.txt
        mencoder mf://@movie.txt -ovc x264 -vf scale=1024:768 -mf type=tiff -o ${PATH_TO_MOVIES}nonaff_vs_affine_z${Z/./_}p${P/./_}e${E/./_}r${R/./_}.avi        
        
        cd ${CWD}

    fi
fi
