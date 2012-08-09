#!/bin/bash
clear
FILENAME="$3"
if [ -z "$FILENAME" ]; then
	FILENAME="movie"
fi
matlab -nojvm -nosplash -r frames\($1,$2\),exit
cd /Applications
./HandBrakeCLI -i $OLDPWD/ycontours_pluslevelset2.avi -o $OLDPWD/"$FILENAME".mp4
cd $OLDPWD
rm ycontours_pluslevelset2.avi
