import os, sys
import numpy as np
import h5py
import json

from seg_watershed import watershed, z_watershed
from seg_util import create_border_mask, writeh5
from waterz import agglomerate

def getScoreFunc(scoreF):
    # aff50_his256
    config = {x[:3]: x[3:] for x in scoreF.split('_')}
    if 'aff' in config:
        if 'his' in config and config['his']!='0':
            return 'OneMinus<HistogramQuantileAffinity<RegionGraphType, %s, ScoreValue, %s>>' % (config['aff'],config['his'])
        else:
            return 'OneMinus<QuantileAffinity<RegionGraphType, '+config['aff']+', ScoreValue>>'
    elif 'max' in config:
        return 'OneMinus<MeanMaxKAffinity<RegionGraphType, '+config['max']+', ScoreValue>>'

def waterz_2d(
            affs,
            useMask=False
            ):
    fragments = z_watershed(np.copy(affs),T_threshes=[300],T_dust=200,T_aff=[0.05,0.8,0.2],T_aff_relative=1-useMask)
    #if useMask:
    #    masks = z_watershed(np.copy(affs),T_threshes=[300],T_dust=200,T_aff_relative=False)
    #    frag_mask = masks > 0
    #    fragments[frag_mask==0]= 0
    return fragments

def waterz(
        affs,
        thresholds,
        output_prefix = './',
        merge_function = None,
        gt = None,
        gt_border = 25/4.0,
        custom_fragments = True,
        discretize_queue = 256,
        fragments_mask = None,
        aff_threshold  = [0.0001,0.9999],
        return_seg = True,
        fragments = None,
        m_history = True,
        useMask = True,
        returnMid = False,
        original_waterz=False):
    #print 'get here.................................................'
    # affs shape: 3*z*y*x
    thresholds = list(thresholds)
    print "waterz at thresholds " + str(thresholds)
    #print affs.shape
    """SRUJAN: original_waterz ultimately refers to a call to mahotas.cwatershed"""
    if original_waterz:
        useMask = False
    if fragments is None:
        if custom_fragments:
            #print np.max(affs[0]) 
            #fragments = z_watershed(np.copy(affs),T_threshes=[150],T_dust=150,T_aff=[0.05,0.8,0.2],T_aff_relative=True)
            if useMask:
                #masks = z_watershed(np.copy(affs),T_threshes=[300],T_dust=200,T_aff_relative=False)  #masks
                fragments = z_watershed(np.copy(affs),T_threshes=[150],T_dust=150,T_aff=[0.05,0.8,0.2],T_aff_relative=False)
            else:
                #print fragments.shape
                if original_waterz:
                    fragments = watershed(affs, 'maxima_distance')
                else:
                    """
                    SRUJAN:
                    This is the case for Haidong's best result. 
                    fragments here refers to watershedded 2D slices of the affinity graph.
                    """
                    fragments = z_watershed(np.copy(affs),T_threshes=[150],T_dust=150,T_aff=[0.05,0.8,0.2],T_aff_relative=True) 
            # for debug
            #print np.max(affs[0])
            if returnMid:
                return fragments  # fragments here means initiate 2d-watershed method
            if fragments_mask is not None:
                fragments[fragments_mask==False] = 0
    #print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    #if useMask:
    #    frag_mask = masks > 0
    #    print 'mask size:-------',np.sum(1-frag_mask)
    #    assert np.sum(frag_mask) > 0 and np.sum(1-frag_mask)>0
    #print fragments.shape
    outs=[]
    if gt is not None and gt_border !=0:
        gt = create_border_mask(gt, gt_border, np.uint64(0))
    #return fragments
    #print fragments
    #print 'agging'
    #for a,b,c in agglomerate(affs,thresholds,gt=gt,aff_threshold_low=aff_threshold[0],aff_threshold_high = aff_threshold[1],return_merge_history = True,fragments=fragments,scoring_function=getScoreFunc(merge_function),discretize_queue=discretize_queue):
    #    print c
    #print '***************************************************************'
    #print 'done'
    """SRUJAN: Haidong did not change the agglomerate function at all"""
    for i,out in enumerate(agglomerate(
            affs,
            thresholds,
            gt = gt,
            aff_threshold_low  = aff_threshold[0],
            aff_threshold_high = aff_threshold[1],
            return_merge_history = True,
            fragments=fragments,
            scoring_function=getScoreFunc(merge_function),
            discretize_queue=discretize_queue)):

        threshold = thresholds[i]
        output_basename = output_prefix+merge_function+'_%.2f'%threshold
        print len(out)
        if gt is not None:
            seg = out[0]
        else:
            seg = out[0]
        # add mask here
        #print np.sum(seg==0)
        #if useMask:
        #    seg[frag_mask==0]= 0
        #print seg.shape
        #print np.sum(seg==0)
        m_history = []
        if return_seg:
            outs.append(seg.copy())
            #m_history.append(out[2])
        else:
            print "Storing segmentation..."
            writeh5(output_basename + '.hdf', 'main', seg)
        if gt is not None:
            metrics = out[1]
            print "Storing record..."
            record = {
                'threshold': threshold,
                'merge_function': merge_function,
                'custom_fragments': custom_fragments,
                'discretize_queue': discretize_queue,
                'voi_split': metrics['V_Info_split'],
                'voi_merge': metrics['V_Info_merge'],
                'rand_split': metrics['V_Rand_split'],
                'rand_merge': metrics['V_Rand_merge'],
            }
            with open(output_basename + '.json', 'w') as f:
                json.dump(record, f)
    if return_seg:
        return outs#,m_history
