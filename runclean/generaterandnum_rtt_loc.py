#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : generaterandnum_rtt_loc.py
#
# Purpose :
#
# Creation Date : 19-01-2019
#
# Last Modified : Fri 19 Apr 2019 05:44:48 PM EDT
#
# Created By : Hongjian Fang: hfang@mit.edu
#
#_._._._._._._._._._._._._._._._._._._._._.*/
def generate_rtt_loc():
	import pandas
	import numpy as np
	import random
	
	#pdcc =pandas.read_csv('dt.cc_syn_loc',header=None,names=['1','2','3','4','5'],delim_whitespace=True)
	pdcc =pandas.read_csv('differential_P.txt',header=None,names=['1','2','3','4','5'],delim_whitespace=True)
	
	nd = len(pdcc)
	nchoose = 100000
	idx = np.sort(random.sample(xrange(nd),nchoose))
	
	df = pdcc.iloc[idx]
	
	with open('tempdata/randomchoosereloc_rtt.dat','wr') as outfile:
		outfile.write(str(nchoose)+'\n')
		outfile.write('\n'.join(
					df['1'].astype(int).astype(str) \
					+' '+df['2'].astype(int).astype(str) \
					+' '+df['3'].astype(int).astype(str) \
					+' '+df['4'].astype(str) \
					+' '+df['5'].astype(str))) 
