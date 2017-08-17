import os
from pickle import load


def pickle_cleaner():
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
		if name.rsplit('.',1)[-1]=='pckl':
            		try:
                    		load(name,'rb')
				print("worked")
            		except:
                		os.remove(name)
				print(os.path.join(root,name))

pickle_cleaner()
