n=$1
cms=$2
./sat_enumeration_brown_dls $n $cms > out_enum_brown_dls_n${n} -cpunum=1
find . -name "brown_dls_n${n}_cover*" -exec cat {} \; > brown_dls_n${n}
rm -f brown_dls_n${n}_cover*
sort -u brown_dls_n${n} > brown_dls_n${n}_nodupl
rm -f brown_dls_n${n}
wc -l brown_dls_n${n}_nodupl
