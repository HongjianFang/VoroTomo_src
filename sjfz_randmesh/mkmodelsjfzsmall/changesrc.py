import numpy as np
dep = 3.0
hori = 0.06
with open('sources.in','r') as src:
	for line in src:
         ls = line.split()
	 if len(ls)>3 and 32<float(ls[1])<35 and 5<float(ls[0])<20:
           print float(ls[0])+dep*np.random.randn(),float(ls[1])+hori*np.random.randn(),float(ls[2])+hori*np.random.randn()
	 else:
	   print line, 
