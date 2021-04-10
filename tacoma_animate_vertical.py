import bpy
import numpy as np
from mathutils import Euler
import math
vertical_armature=bpy.data.objects['Armature.001']
A=np.loadtxt('/home/kmd/Dropbox/research/Dynamics/bridge/video_assets/txtdata/tacoma_sync.mcol')
tower_z0=vertical_armature.pose.bones['TowerNorthAnim'].location[1]
h0=abs(vertical_armature.pose.bones['TowerNorthAnim'].head[1]-vertical_armature.pose.bones['TowerSouthAnim'].head[1])
deck_z0=vertical_armature.pose.bones['Deck'].location[0]
for i,col in enumerate(A[:,:]):
    t=col[0]
    vertical_armature.pose.bones['TowerSouthAnim'].location[1]=tower_z0+col[1]/20.0
    vertical_armature.pose.bones['TowerSouthAnim'].keyframe_insert(data_path='location', frame=i)
    vertical_armature.pose.bones['TowerNorthAnim'].location[1]=tower_z0+col[2]/20.0
    vertical_armature.pose.bones['TowerNorthAnim'].keyframe_insert(data_path='location', frame=i)
    length=(col[2]-col[1])
    #vertical_armature.pose.bones['Deck'].rotation_euler=Euler([0,0,math.asin(length/h0)])
    vertical_armature.pose.bones['Deck'].location[0]=deck_z0+col[3]/400.0
  
    #vertical_armature.pose.bones['Deck'].keyframe_insert(data_path='rotation_euler', frame=i)
    vertical_armature.pose.bones['Deck'].keyframe_insert(data_path='location', frame=i)
    