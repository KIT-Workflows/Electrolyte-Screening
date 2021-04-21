#!/bin/bash

module purge
module load turbomole/7.4.1

{%- set general = wano["General"] -%}
{%- set struct = general["Molecular structure"] -%}
{%- set eht = general["Starting orbitals (Extended Hueckel Guess)"] -%}
{%- set method = general["Method"] %}

{% if struct["Structure file type"] == "Turbomole coord" -%}
mv structure_file coord_0
{% else -%}
x2t structure_file > coord_0
{% endif %}

define << EOF > define.out



a coord_0
{% if struct["Use symmetry"] == True -%}
desy {{ struct["Threshold"] }} 
{% endif %}
{% if struct["Use internal coordinates"] == True -%}
ired
*
{% else %}
*
no
{% endif %}
b all {{ general["Basis set"]["Basis set"] }}
*
eht

{{ eht["Charge"] | int }}
{% if eht["Multiplicity"] > 2 -%}
n
u {{ eht["Multiplicity"]-1 }}
*
{% else -%}

{% endif %}


{% if method["Method"] == "DFT" -%}
ri
on
m 1000
*
dft
on
func {{ method["Functional"] }}
grid m4
*
{% if method["Dispersion correction"] != "None" -%}
dsp
{% if method["Dispersion correction"] == "D2" -%}
old
{% elif method["Dispersion correction"] == "D3" -%}
on
{% else -%}
bj
{% endif %}
*
{% endif %}
{% else %}
rijk
on
m 1000
*
{% endif %}
{% if method["Method"] == "Post-HF" -%}
cc
cbas
*
{% if method["Frozen Core"] != "None" -%}
freeze
{% if method["Frozen Core"] == "Manual" -%}
fp -1000.0
f 1-{{ method["Number of frozen orbitals"] }}
{% endif %}
*
{% endif %}
memory
1000
ricc2
model {{ method["Correlation treatment"] }}
*
*
{% endif %}
*
EOF

ridft > ridft.out
eiger > eiger.out
{% if method["Method"] == "Post-HF" -%}
ccsdf12 > ccsdf12.out
{%- endif %}

exit 0
