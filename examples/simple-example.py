#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import cv2
import pydegensac
from time import time
from copy import deepcopy

#Now helper function for running homography RANSAC
def verify_cv2(kps1, kps2, tentatives, th = 4.0 , n_iter = 2000):
    src_pts = np.float32([ kps1[m.queryIdx].pt for m in tentatives ]).reshape(-1,2)
    dst_pts = np.float32([ kps2[m.trainIdx].pt for m in tentatives ]).reshape(-1,2)
    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, th, 0.99, n_iter)
    print ('cv2 found {} inliers'.format(int(deepcopy(mask).astype(np.float32).sum())))
    return H, mask

def verify_pydegensac(kps1, kps2, tentatives, th = 4.0,  n_iter = 2000):
    src_pts = np.float32([ kps1[m.queryIdx].pt for m in tentatives ]).reshape(-1,2)
    dst_pts = np.float32([ kps2[m.trainIdx].pt for m in tentatives ]).reshape(-1,2)
    H, mask = pydegensac.findHomography(src_pts, dst_pts, th, 0.99, n_iter)
    print ('pydegensac found {} inliers'.format(int(deepcopy(mask).astype(np.float32).sum())))
    return H, mask

def verify_cv2_fundam(kps1, kps2, tentatives, th = 1.0 , n_iter = 10000):
    src_pts = np.float32([ kps1[m.queryIdx].pt for m in tentatives ]).reshape(-1,2)
    dst_pts = np.float32([ kps2[m.trainIdx].pt for m in tentatives ]).reshape(-1,2)
    F, mask = cv2.findFundamentalMat(src_pts, dst_pts, cv2.RANSAC, th, 0.999, n_iter)
    print ('cv2 found {} inliers'.format(int(deepcopy(mask).astype(np.float32).sum())))
    return F, mask

def verify_pydegensac_fundam(kps1, kps2, tentatives, th = 1.0,  n_iter = 10000):
    src_pts = np.float32([ kps1[m.queryIdx].pt for m in tentatives ]).reshape(-1,2)
    dst_pts = np.float32([ kps2[m.trainIdx].pt for m in tentatives ]).reshape(-1,2)
    F, mask = pydegensac.findFundamentalMatrix(src_pts, dst_pts, th, 0.999, n_iter, enable_degeneracy_check= True)
    print ('pydegensac found {} inliers'.format(int(deepcopy(mask).astype(np.float32).sum())))
    return F, mask

if __name__ == '__main__':
    img1 = cv2.cvtColor(cv2.imread('img/v_dogman/1.ppm'), cv2.COLOR_BGR2RGB)
    img2 = cv2.cvtColor(cv2.imread('img/v_dogman/6.ppm'), cv2.COLOR_BGR2RGB)
    # SIFT is not available by pip install, so lets use AKAZE features
    det = cv2.AKAZE_create(descriptor_type = 3, threshold=0.00001)
    kps1, descs1 = det.detectAndCompute(img1,None)
    kps2, descs2 = det.detectAndCompute(img2,None)
    bf = cv2.BFMatcher()
    matches = bf.knnMatch(descs1,descs2, k=2)
    matchesMask = [False for i in range(len(matches))]
    # SNN ratio test
    for i,(m,n) in enumerate(matches):
        if m.distance < 0.9*n.distance:
            matchesMask[i]=True
    tentatives = [m[0] for i, m in enumerate(matches) if matchesMask[i] ]

    th = 4.0
    n_iter = 2000
    t=time()
    print ("Running homography estimation")
    cv2_H, cv2_mask = verify_cv2(kps1,kps2,tentatives, th, n_iter )
    print ("OpenCV runtime {0:.5f}".format(time()-t), ' sec')
    t=time()
    cmp_H, cmp_mask = verify_pydegensac(kps1,kps2,tentatives, th, n_iter)
    print ("pydegensac runtime {0:.5f}".format(time()-t), ' sec')
    print ("H = ", cmp_H)
    th = 0.5
    n_iter = 50000
    print ("Running fundamental matrix estimation")

    t=time()
    cv2_H, cv2_mask = verify_cv2_fundam(kps1,kps2,tentatives, th, n_iter )
    print ("OpenCV runtime {0:.5f}".format(time()-t), ' sec')

    t=time()
    cmp_H, cmp_mask = verify_pydegensac_fundam(kps1,kps2,tentatives, th, n_iter)
    print ("pydegensac {0:.5f}".format(time()-t), ' sec')
    print ("F = ", cmp_H)
