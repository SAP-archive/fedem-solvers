#!/bin/sh
# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Bash script that post-processes the two frs-files resulting from running
# the respons_pos program, into a single file response_pos.frs

echo '' | tee -a response_pos_1.frs >> gage_pos_2.frs
sed '1,/Physical time/d;/^\[/,$ d' gage_pos_2.frs > gage_1.frs
grep -a '^\[' gage_pos_2.frs > gage_2.frs
grep -a '^{' gage_pos_2.frs > gage_3.frs
egrep -a -n '^\[1;|^DATABLOCKS:|^DATA:' response_pos_1.frs |\
awk -F: 'BEGIN{i=1;}{print$1-1" r gage_"i++".frs"}' > insert.sed
sed -f insert.sed response_pos_1.frs |\
sed '/^{/s/ *[0-9]; *[0-9];""/;;/g' > response_pos.frs
rm response_pos_1.frs gage_*.frs insert.sed
