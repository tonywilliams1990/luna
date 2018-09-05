#!/bin/sh
# Script for commiting Luna to GitHub
echo "Commit script"

make clean

# Move unnecessary files
#mv ./DATA ~/Desktop/DATA_temp
#mkdir ./DATA
#cp ~/Desktop/DATA_temp/.gitkeep ./DATA/.gitkeep
#rm -rf .sconf_temp

# Setup commit message
if [ $# -eq 0 ] ;
then
  message="empty"
else
  message=$1
fi

echo $message

# GitHub stuff
git add .
git commit -m $message
git push origin master

# Move the files back
#rm -r DATA
#mv ~/Desktop/DATA_temp ./DATA

exit 0
