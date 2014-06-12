#!/bin/bash

mysql -uhuebeli_vmc -pkerberos -hhuebeli.mysql.db.hostpoint.ch huebeli_vmc -ss -e "$1" > $2.out

zip -j $2.zip $2.out

