#!/bin/bash

outfile=$2.out

mysql -uhuebeli_vmc -pkerberos -hhuebeli.mysql.db.hostpoint.ch huebeli_vmc -ss -e "$1" > $outfile

sed s/NULL/0/g $outfile > ${outfile}0
mv ${outfile}0 $outfile

zip -j $2.zip $2.out

