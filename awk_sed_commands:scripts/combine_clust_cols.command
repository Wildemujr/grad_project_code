#!/bin/bash

fh1=$1
fh2=$2
fh3=$3

read -p $'\nWhat is the name of the output file?\n' f_na

gawk '{ if (NR!=1) {print $1" "$2} }' ${1}.txt > ${f_na}.1
gawk '{ if (NR!=1) {print $2} }' ${2}.txt > ${f_na}.2
gawk '{ if (NR!=1) {print $2} }' ${3}.txt > ${f_na}.3
paste -d' ' ${f_na}.1 ${f_na}.2 > ${f_na}.12.txt
paste -d' ' ${f_na}.12.txt ${f_na}.3  > ${f_na}.123.txt
awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= 3; print sum}' ${f_na}.123.txt > ${f_na}.avg
paste -d' ' ${f_na}.123.txt ${f_na}.avg > ${f_na}.fin.txt
sed '1i Frame_Number Distance_1 Distance_2 Distance_3 Average_Distance' ${f_na}.fin.txt > ${f_na}.fin.fin.txt
sed -E 's/[[:space:]]+/,/g' ${f_na}.fin.fin.txt > ${f_na}.avg.csv
# gawk -F ',' '{ print $1,$2,$3,$4 }' ${f_na}.123 > ${f_na}.avg.csv
rm -f ${f_na}.1 ${f_na}.2 ${f_na}.3 ${f_na}.12.txt ${f_na}.123.txt \
	${f_na}.fin.txt ${f_na}.fin.fin.txt ${f_na}.avg
 