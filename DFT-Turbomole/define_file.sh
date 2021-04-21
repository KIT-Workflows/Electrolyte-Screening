#!/bin/bash

module purge
module load turbomole/7.4.1

x2t structure_file > coord

define << EOF > define.out



a coord
*


ired
*

b all def-SV(P)
*
eht

0


dft
on
func pbe
grid m4
dsp
old
*

ri
on
m 2000
*


EOF

ridft > ridft.out
eiger > eiger.out

jobex -dscf -c 500 > jobex.out

exit 0
