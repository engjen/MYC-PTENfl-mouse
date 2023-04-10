# 2020-05-05 Need to threshold with Napari!
# author engje
# ipython --gui=qt
#%run 20200504_JPTMAs_napari.py
import napari
import os
import skimage
from skimage import io, measure
import numpy as np
import copy
import pandas as pd
import tifffile

#paths
codedir = 'C:\\Users\\engje\\Documents\\Data\\cmIF_2021-07-07_RS-mTMA'

regdir = f'{codedir}\\Cropped'

os.chdir('..')
from mplex_image import visualize as viz
from mplex_image import process
from mplex_image import analyze

def add_points(s_crop,df,s_slide,i_scale = 2):
    x_corr = int(s_crop.split('x')[1].split('y')[0])
    y_corr = int(s_crop.split('x')[1].split('y')[1])
    x_arr = (df.loc[df.index.str.contains(s_slide),['DAPI_X']] - x_corr)/i_scale
    y_arr = (df.loc[df.index.str.contains(s_slide),['DAPI_Y']] - y_corr)/i_scale
    x_arr['DAPI_Y'] = y_arr
    points = np.array(x_arr.loc[:,['DAPI_Y','DAPI_X']])
    return(points,x_arr)

#load positive and intensity data
os.chdir(codedir)

s_out = '20211111_RS-mTMA-5'
df_man = pd.read_csv(f'{s_out}_GatedPositiveCellNames.csv',index_col=0)
df_pos = analyze.celltype_to_bool(df_man,'celltype')
df_pos['slide_scene'] = [item.split('_cell')[0] for item in df_pos.index]

df_mi = pd.read_csv(f'20211109_RS-mTMA-5_FilteredMeanIntensity.csv',index_col=0)

# weird ones
s_slide = 'RS-mTMA-5_sceneA05' #high foxp3 bg, some float
s_slide = 'RS-mTMA-5_sceneE03' #corner extra
s_slide = 'RS-mTMA-5_sceneE07' #folded tissue
s_slide = 'RS-mTMA-5_sceneH06' #bright floating tissue
s_slide = 'RS-mTMA-5_sceneI01' #bright floating tissue
s_slide = 'RS-mTMA-5_sceneI11' #epcam, ecad, CK5 neg (no tumor), foxp3 floating tissue
s_slide = 'RS-mTMA-5_sceneH11' # circel imaging artifact (pmyc)
s_slide = 'RS-mTMA-5_sceneI07' #dried out bottom? pmyc weird
s_slide = 'RS-mTMA-5_sceneI02' #cd31 bg, epcam bright on top, dim on bottom. dried out bottom?

#new outliers
s_slide = 'RS-mTMA-5_sceneA05' #lots of CD45, Ecad blurry/neg, CK5 and Epcam negative, sox9 pos
s_slide = 'RS-mTMA-5_sceneE09' #cd45 artifact
s_slide = 'RS-mTMA-5_sceneF03' #exclude upper right corner
s_slide = 'RS-mTMA-5_sceneF04' #axclude scattered cells on edge
#s_slide = 'RS-mTMA-5_sceneG09' #circle artifact, AF area
s_slide = 'RS-mTMA-5_sceneH06' #pmyc weird, not sure where

#pretty
s_slide = 'RS-mTMA-5_sceneA09'
#s_slide = 'RS-mTMA-5_sceneB09'
#s_slide = 'RS-mTMA-5_sceneB07'

#stroma poor
#s_slide = 'RS-mTMA-5_sceneI08'  messed up with clustering
s_slide = 'RS-mTMA-5_sceneA01' #true - some stroma; wrong annotation by xiaoyan?
s_slide = 'RS-mTMA-5_sceneA04' #true - lots of stroma wrong annot?

#s_slide = 'RS-mTMA-5_sceneG05' #tru - has lots of stroma - wrong annot?

#cluden low?
s_slide = 'RS-mTMA-5_sceneI11' #no tumor; #Ecad-, CK5-, EpCAM-; (claudin low)? 
#s_slide = 'RS-mTMA-5_sceneG09' #Ecad-, CK5-, EpCAM-; (claudin low)? 
s_slide = 'RS-mTMA-5_sceneH11' #also ecad neg, epcam neg (claudin low)? 
#s_slide = 'RS-mTMA-5_sceneA05' #could drop - Ecad stainifn weird
#s_slide = 'RS-mTMA-5_sceneF03' #very weird : pMYC negative and claugin low; drop.
#s_slide = 'RS-mTMA-5_sceneA09'
#s_slide = 'RS-mTMA-5_sceneB07'
#s_slide = 'RS-mTMA-5_sceneC08'
#s_slide = 'RS-mTMA-5_sceneF02'
#s_slide = 'RS-mTMA-5_sceneC01' #good IR
#s_slide = 'RS-mTMA-5_sceneG06' #good SP
s_slide = 'RS-mTMA-5_sceneD07' #btter SP
#s_slide = 'RS-mTMA-5_sceneI10'
#s_slide = 'RS-mTMA-5_sceneH03'
s_slide = 'RS-mTMA-5_sceneF11'

#check

#load images
print(s_slide)
os.chdir(regdir)
for s_file in os.listdir():
    if s_file.find(s_slide) > -1:
        s_crop = s_file.split('_')[2]
viewer = napari.Viewer()
label_image = viz.load_crops(viewer,s_crop,s_slide)

#show positive dots
ls_cell = ['tumor','stromal','immune']
df_scene = df_mi[(df_mi.slide_scene==s_slide)]
d_cell = {'tumor':'blue','stromal':'orange','immune':'green'}

for s_cell in ls_cell: 
    ls_index = df_pos[df_pos.loc[:,s_cell]].index
    points,x_arr = add_points(s_crop,df_scene[df_scene.index.isin(ls_index)],s_slide)
    viewer.add_points(points,name=s_cell,face_color=d_cell[s_cell])


# centroids
'''
points,x_arr = add_points(s_crop,df_mi,s_slide)
viewer.add_points(points)
'''
#draw shape
'''
verts = viewer.layers['Shapes'].data[0]
b_poly = measure.points_in_poly(points, verts)
point_properties = {
    'label': df_mi.loc[df_mi.index.str.contains(s_slide),['cellid']],
    'in_poly' : b_poly
}
points_layer = viewer.add_points(points, properties=point_properties, face_color='in_poly',size=10)
'''
#select mask
'''
#inpoly
df_exclude = x_arr.loc[b_poly]

#out of poly
df_exclude = x_arr.loc[~b_poly]

df_exclude = df_exclude.append(x_arr.loc[b_poly])
df_exclude = df_exclude.append(x_arr.loc[~b_poly])

df_exclude.to_csv(f'exclude_{s_slide}.csv')
'''


#show positive results
'''
ls_cell = ['tumor','stromal']
df_scene = df_pos[df_pos.slide_scene==s_slide]
#df_scene = df_pos[df_pos.slide_scene==f"{s_slide.split('-Scene-')[0]}_scene{s_slide.split('-Scene-')[1]}"]

for s_cell in ls_cell: 
    label_image_cell = viz.pos_label(viewer,df_scene,label_image,s_cell)
'''
os.chdir(codedir)


