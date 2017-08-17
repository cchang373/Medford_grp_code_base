import os
from tempfile import mkstemp
from shutil import move
from os import remove, close

#from subprocess import call
cur_dir=os.path.relpath(".","..")
import os
import sys
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
sys.path.insert(0,"/nv/hp13/bcomer3/shared/espresso_rutile/tools")
from Change_Cores import custom_replace

textToSearch = "joe-6-ge"
textToReplace = "joe"


for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if name == "run.sh"or name =='run_vib.py':
            print(os.path.join(root, name))
            custom_replace(os.path.join(root, name),textToSearch, textToReplace)
        else:
            pass


