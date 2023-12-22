#!/bin/sh
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# This script builds the html-documentation for fedempy using sphinx.
# It will do the same as "make.bat html" in a DOS shell, but will
# in addition insert the correct version tags in the HOWTO install page.

if command -v sphinx-build &> /dev/null; then
  cd `dirname $0`
  sphinx-build -M html source build
else
  echo The 'sphinx-build' command was not found. Make sure you have Sphinx installed.
  exit 1
fi

ftpyver=`cat ../version.txt`
buildno=`sed "s/^.*\.//" ../version.txt`
version=`sed "1 s/\"//g;1 s/.*$/& (build $buildno)/" ../../fedem-foundation/src/Admin/version.h`
sed -i "s/#FEDEMPY_VERSION#/$ftpyver/;s/#FEDEM_VERSION#/$version/" build/html/howto-install.html
