#!/bin/sh
# Helper script to insert the current Fedem version with build number into the generated HOWTO documentation.
# $1 = the html-file to edit
# $2 = the file holding the current Fedem version
# $3 = the file holding the current build number (using the last digit only)
if [ $# -lt 3 ]; then exec echo "usage: $0 <htmlfile> <versionfile> <buildfile>"; fi
ftpyver=`cat $3`
buildno=`sed "s/^.*\.//" $3`
version=`sed "1 s/\"//g;1 s/.*$/& (build $buildno)/" $2`
sed -i "s/#FEDEMPY_VERSION#/$ftpyver/;s/#FEDEM_VERSION#/$version/" $1
