# -*- coding: utf-8 -*-
"""
Created on Mon May 29 11:03:34 2017

@author: benjamin
"""
from ase import io
import os
from mapping import mapping_dict
kwargs = {
  #  'rotation'      : rot, # text string with rotation (default='' )
    'radii'         : .85, # float, or a list with one float per atom
    'colors'        : None,# List: one (r, g, b) tuple per atom
    'show_unit_cell': 0,   # 0, 1, or 2 to not show, show, and show all of cell
    }

# Extra kwargs only available for povray (All units in angstrom)
kwargs.update({
    'run_povray'   : True, # Run povray or just write .pov + .ini files
    'display'      : False,# Display while rendering
    'pause'        : True, # Pause when done rendering (only if display)
    'transparent'  : False,# Transparent background
    'canvas_width' : None, # Width of canvas in pixels
    'canvas_height': 750, # Height of canvas in pixels 
    'camera_dist'  : 50.,  # Distance from camera to front atom
    'image_plane'  : None, # Distance from front atom to image plane
    'camera_type'  : 'perspective', # perspective, ultra_wide_angle
    'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
    'area_light'   : [(2., 3., 40.), # location
                      'White',       # color
                      .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
    'background'   : 'White',        # color
    'textures'     : None, # Length of atoms list of texture names
    'celllinewidth': 0.1,  # Radius of the cylinders representing the cell
    })
def make_pics(path,fancy=True):
    name = path.split(os.sep)[-3]+'-'+path.split(os.sep)[-2]+'-'+path.split(os.sep)[-4]
    if name in mapping_dict.keys():
        slb = io.read(path+'/converged_slab.traj')
        if fancy == True:
            io.write('/home/benjamin/AJ_Compounds/'+name+'.pov',slb,rotation='-45x',**kwargs)
            #io.write('/home/benjamin/paper1/N2_fixation_paper/figures/compounds/'+name+'-xside.pov',slb,rotation='-90x',**kwargs)
            #io.write('yside.png',slb,rotation='-90y')
            #image = cv2.imread('yside.png')
            print('lol')
        else:
            io.write('/home/benjamin/paper1/N2_fixation_paper/figures/compounds/'+name+'-top.png',slb)
            io.write('/home/benjamin/paper1/N2_fixation_paper/figures/compounds/'+name+'-xside.png',slb,rotation='-90x')

    
for root, dirs, files in os.walk(".", topdown=False):
    for name in dirs:
        #print(os.path.join(root,name))
        esp_pres = [x for x in root.split(os.sep) if x=="esp.log"]
        non_vib = esp_pres==["esp.log"]
        slab_conv_pres = [x for x in root.split(os.sep) if x=="slab_convergence"]
        non_save = [x for x in root.split(os.sep) if x=="calc.save"] == []
        non_slab_conv = slab_conv_pres ==[]
        #print os.path.isfile(name+"/converged_slab.traj") and non_vib==False and non_slab_conv == True
        if os.path.isfile(os.path.join(root,name)+"/converged_slab.traj") and non_vib==False and non_slab_conv == True:
            #print name
            make_pics(os.path.join(root,name))