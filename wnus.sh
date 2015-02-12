nline=`wc -l nu_vol.out | tr -d '[:alpha:][= =][=.=][=_=]'`
nline=`echo "$nline/2" | bc`
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' nu_vol.out`
echo "Half Nu_conv: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nu_vol.out`
echo "Full Nu_conv: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nu_plate.out`
echo "Full Nu_bot: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' nu_plate.out`
echo "Full Nu_top: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' nu_diss.out`
echo "Full Nu_epsu: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' nu_diss.out`
echo "Full Nu_epst: $nu"

