import sys
import os
import getopt
from random import randint

if __name__=='__main__':
	for i in range(19, 25):
		os.system('nohup ./our/CVSPsearch Vmaps/map{} Pmaps/phy -l 50 -f -c 0 500 -n 100 -o Outputs/Exp1-{} &> Outputs/Exp1-{}-nohup &'.format(i,i,i))