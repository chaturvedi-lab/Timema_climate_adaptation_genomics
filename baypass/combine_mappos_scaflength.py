#!/usr/bin/python2.7
import csv
import sys

## ------------------------------------------------- ##
## Author: Samridhi Chaturvedi
## Copyright: Copyright 2019
## Credits: [""]
## License: GPL3
## Version: 0.0.1
## Maintainer: Samridhi Chaturvedi
## Email: samridhi.chaturvedi@gmail.com
## -------------------------------------------- ##

'''
Usage:
        python2.7 combine_mappos_scaflength.py scaffold_length_file.txt mappos_file.txt outfile.txt
ot.txt 

Prerequisite:
Make sure the analyses file is space separated.

'''

##read in the files
scaflenfile = open(sys.argv[1],'r')
mapposfile = open(sys.argv[2], 'r')
outfile = open(sys.argv[3],'w')
outheader_str = "scaffold\tposition\tlength\n"
outfile.write(outheader_str)
                        
scafpos = dict()
next(scaflenfile)
for line in scaflenfile:
        line = line.strip('\n')
        line = line.split(' ')
        #print line
        keyannot = line[0]
        #print keyannot
        scafpos[keyannot] = line
        #print scafpos

for aline in mapposfile:
        aline = aline.strip('\n')
        aline = aline.split('\t')
        #print aline
        keymp = aline[0]
        w = aline
        w.append('Missing')
        if keymp in scafpos:
                match = scafpos[keymp]
                #print match
                w[2] = match[1]
                #print w
                out = '\t'.join(w)
                out_str = str(out)
                #print out_str
                outfile.write(out_str + '\n')
