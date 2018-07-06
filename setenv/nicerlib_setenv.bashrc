#!/bin/bash 

NAME=$(scutil --get ComputerName)
if [ $NAME = 'vicuna' ]; then
	echo '...setting for machine of "vicuna"'
	export NICER_DATA_PATH=/Users/enoto/work/nicer/data
	export NICER_PUBLIC_DATA_PATH=/Users/enoto/work/nicer/data/obs	
	export NICER_SOFT_PATH=/Users/enoto/work/soft/nicerlib/nicerlib
	export NICER_RESP_PATH=/Users/enoto/work/niresp
elif [ $NAME = 'nebula' ]; then
	echo '...setting for machine of "nebula"'
	export NICER_DATA_PATH=/Users/enoto/work/drbv2/nicer/data
	export NICER_PUBLIC_DATA_PATH=$NICER_DATA_PATH/obs
	export NICER_SOFT_PATH=/Users/enoto/work/drbv2/nicer/soft/nicerlib
	export NICER_RESP_PATH=/Users/enoto/work/niresp
elif [ $NAME = 'llama' ]; then
	echo '...setting for machine of "fireant"'
	export NICER_DATA_PATH=/Users/enoto/work/nicer/data
	export NICER_PUBLIC_DATA_PATH=/Users/enoto/work/nicer/data/obs
	export NICER_SOFT_PATH=/Users/enoto/work/nicer/soft/nicerlib
	export NICER_RESP_PATH=/Users/enoto/work/niresp
else
	echo 'no corresponding computer setup.'
fi

export NICER_BGDMODEL_PATH=/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd

export PYTHONPATH=$NICER_SOFT_PATH:$PYTHONPATH

export PATH=$NICER_SOFT_PATH/scripts:$PATH
export PATH="/Users/enoto/Dropbox/enoto/research/nicer/nisoft/ver3.00/local":$PATH
export PATH=$NICER_SOFT_PATH/scripts/repository:$PATH
export PATH=$NICER_SOFT_PATH/scripts/mitbgd/3C50:$PATH
export PATH=$NICER_SOFT_PATH/scripts/process/proc_kyoto180608:$PATH
export PATH=$NICER_SOFT_PATH/scripts/speccorr:$PATH
export PATH=$NICER_SOFT_PATH/scripts/misc:$PATH
export PATH=$NICER_SOFT_PATH/scripts/ql:$PATH

echo ComputerName    = $NAME
echo NICER_DATA_PATH = $NICER_DATA_PATH
echo NICER_PUBLIC_DATA_PATH = $NICER_PUBLIC_DATA_PATH
echo NICER_SOFT_PATH = $NICER_SOFT_PATH
echo NICER_RESP_PATH = $NICER_RESP_PATH
echo NICER_BGDMODEL_PATH = $NICER_BGDMODEL_PATH
echo PYTHONPATH = $PYTHONPATH
