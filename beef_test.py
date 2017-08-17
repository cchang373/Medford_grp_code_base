from Beef_data import beef
import numpy as np
l = beef([1.1,np.array([1,1])])
k = beef([2.1,np.array([2,2])])
j = beef([1.1,np.array([1,1])])
p = l+k-j
l+=j
k-=p
l = l+1
l[0]
l[:-1]
l[-1]
l.error_bar()
print(l)
h = l
h
h = l.c()
h
h.add_path('lol')
n = h**5
n.s()

