#!/bin/sh
for i in `ls -1 f*tex`; do
    rubber $i
    dvisvgm --no-fonts --scale=1.5 `basename $i tex`dvi
    rubber --clean $i
done

#for i in `ls -1 *pdf`; do
#    convert -density 150 $i `basename $i pdf`png
#done


