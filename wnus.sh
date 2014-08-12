nline=`wc -l nusse.out | tr -d '[:alpha:][= =][=.=]'`
nline=`echo "$nline/2" | bc`
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' nusse.out`
echo "Half Nu_conv: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nusse.out`
echo "Full Nu_conv: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nusse2.out`
echo "Full Nu_bot: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' nusse2.out`
echo "Full Nu_top: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nusse3.out`
echo "Full Nu_epsu: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' nusse3.out`
echo "Full Nu_epst: $nu"

