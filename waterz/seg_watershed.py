import numpy as np
import mahotas
from scipy import ndimage
from zwatershed import zwatershed

def get_seeds(boundary, method='grid', next_id = 1,
             seed_distance = 10):
    if method == 'grid':
        height = boundary.shape[0]
        width  = boundary.shape[1]

        seed_positions = np.ogrid[0:height:seed_distance, 0:width:seed_distance]
        num_seeds_y = seed_positions[0].size
        num_seeds_x = seed_positions[1].size
        num_seeds = num_seeds_x*num_seeds_y
        seeds = np.zeros_like(boundary).astype(np.int32)
        seeds[seed_positions] = np.arange(next_id, next_id + num_seeds).reshape((num_seeds_y,num_seeds_x))

    if method == 'minima':
        minima = mahotas.regmin(boundary)
        seeds, num_seeds = mahotas.label(minima)
        seeds += next_id
        seeds[seeds==next_id] = 0

    if method == 'maxima_distance':
        distance = mahotas.distance(boundary<0.5)
        maxima = mahotas.regmax(distance)
        seeds, num_seeds = mahotas.label(maxima)
        seeds += next_id
        seeds[seeds==next_id] = 0
    return seeds, num_seeds

def watershed(affs, seed_method, use_mahotas_watershed = True):
    affs_xy = 1.0 - 0.5*(affs[1] + affs[2])
    depth  = affs_xy.shape[0]
    fragments = np.zeros_like(affs[0]).astype(np.uint64)
    next_id = 1
    for z in range(depth):
        seeds, num_seeds = get_seeds(affs_xy[z], next_id=next_id, method=seed_method)
        if use_mahotas_watershed:
            fragments[z] = mahotas.cwatershed(affs_xy[z], seeds)
        else:
            fragments[z] = ndimage.watershed_ift((255.0*affs_xy[z]).astype(np.uint8), seeds)
        next_id += num_seeds

    return fragments

def z_watershed(affs, T_threshes=[300], T_aff=[0.05,0.8,0.2], T_aff_relative=True, T_dust=150, T_merge=0):
    print 'using 2-d zwatershed on affs:',affs.shape
    _,depth,width,height = affs.shape
    fragments = np.zeros_like(affs[0]).astype(np.uint64)
    #fragments = np.array(zwatershed(affs,T_threshes,T_aff,T_aff_relative,T_dust,T_merge)[0][0])
    #print fragments.shape
    for z in range(depth):
        """SRUJAN: generate z-slices of the affinity graph. fraction is one slice"""
        fraction = np.reshape(affs[:,z,:,:],(3,1,width,height))
        """SRUJAN: This essentially sets all affinities in z-direction to zero"""
        #fraction[0] = 0
        """SRUJAN: zwatershed is an .so object generated from zwatershed.cpp in the source"""
        output = np.array(zwatershed(fraction,T_threshes,T_aff,T_aff_relative,T_dust,T_merge)[0][0])
        """SRUJAN: add an ofset to all nonzero labels in this slice
        so that segment IDs do not repeat from slice to slice"""
        out_mask = output>0
        output += np.max(fragments) * out_mask  # need to re-id here
        #output += np.max(fragments)+1
        fragments[z] = np.reshape(output,(width,height))
    return fragments
