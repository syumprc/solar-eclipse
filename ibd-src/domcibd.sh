#!/bin/sh
rm -f ibdmc.err
if [ "$SOLAR_BIN" = "" ]
then
echo Environment variable SOLAR_BIN not set in domcibd. > ibdmc.err
exit 1 
fi
if [ $# -ne 4 ] ; then
echo "usage: domcibd ibdfile #reps max1 showstat"
exit 1
fi
ibdfil=$1
nreps=$2
max1=$3
showstat=$4
if [ $showstat = 'y' ] ; then
if [ ! -f /usr/bin/printf ] ; then
showstat=n
fi
fi
if [ ! -f translat.tab ] ; then
echo MENDEL input file \"translat.tab\" not found. > ibdmc.err
exit 1
fi
if [ ! -f ibd.loc ] ; then
echo MENDEL input file \"ibd.loc\" not found. > ibdmc.err
exit 1
fi
if [ ! -f ibd.bat ] ; then
echo MENDEL input file \"ibd.bat\" not found. > ibdmc.err
exit 1
fi
rm -f $ibdfil
nfam=`grep FAM translat.tab | wc -l`
i=0
while [ $i -lt $nfam ]
do
i=`expr $i + 1`
fam=`echo $i | awk '{printf "%05d", $0}'`
if [ $showstat = 'y' ] ; then
printf "pedigree%5d " $i >/dev/tty
fi
awk 'BEGIN{prt=1}
	/FAM/{prt=0}
	/FAM'$fam'/{prt=1}
	{if(prt)print $0}' translat.tab > ibd.in
$SOLAR_BIN/ibdmc $nreps $max1 > ibdmc.err 2>/dev/null
if [ $? != 0 ] ; then
if [ ! -s ibdmc.err ] ; then
echo ibdmc: Not enough memory. > ibdmc.err
fi
exit 1
fi
sed 1d ibd.mat >> $ibdfil
rm -f risk.mem lnlik.mem
if [ $showstat = 'y' ] ; then
printf "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" >/dev/tty
fi
done
gzip -f $ibdfil
rm -f ibdmc.err ibd.mat ibd.in ibd.out ibd.ped ibd.rsk
