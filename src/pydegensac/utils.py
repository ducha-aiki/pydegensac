import pydegensac

import numpy as np
import math

try:
    import cv2
    OPENCV_HERE = True
except:
    OPENCV_HERE = False

import warnings


error_type_dict_homography = {"sampson": 0,
                   "symm_sq_max": 1,
                   "symm_max": 2,
                   "symm_sq_sum": 3,
                   "symm_sum": 4}

error_type_dict_fundamental = {"sampson": 0,
                   "symm_epipolar": 1}

def convert_cv2_kpts_to_xyA(kps):
    """
    Converts OpenCV keypoints into transformation matrix
    and pyramid index to extract from for the patch extraction 
    """
    num = len(kps)
    out = np.zeros((num,6)).astype(np.float64)
    for i, kp in enumerate(kps):
        out[i,:2] = kp.pt
        s = kp.size
        a = kp.angle
        cos = math.cos(a * math.pi / 180.0)
        sin = math.sin(a * math.pi / 180.0)
        out[i,2] = s*cos
        out[i,3] = s*sin
        out[i,4] = -s*sin
        out[i,5] = s*cos
    return out

def convert_and_check(kps1):
    if type(kps1) is np.ndarray:
        sh = kps1.shape
        err_message = ValueError("Keypoints should be list of cv2.KeyPoint \
                             or numpy.array [Nx2] or [Nx6]. N>=4  \
                             Got shape of {} with shape instead".format(str(sh))) 
        if len(sh) != 2:
            raise err_message
        num, dim = sh
        if (dim != 2) and (dim!=6):
            raise err_message
        if (num < 4):
            raise err_message 
        out = kps1.astype(np.float64)
    elif type(kps1) is list:
        if OPENCV_HERE:
            if type(kps1[0]) is not cv2.KeyPoint:
                raise ValueError("Keypoints should be list of cv2.KeyPoint \
                                or numpy.array [Nx2] or [Nx6]. N>=4  \
                                Got input of list of type {}".format(str(type(kps1[0])))) 
            else:
                out = convert_cv2_kpts_to_xyA(kps1)
        else:
            raise ValueError("Cannot import cv2. Please, install or pass np.arrays instead")
    else:
        raise ValueError("Keypoints should be list of cv2.KeyPoint \
                             or numpy.array [Nx2] or [Nx6]. N>=4  \
                             Got input of type {}".format(str(type(kps1)))) 
    return out


def findHomography(pts1_,
                   pts2_, 
                   px_th = 1.0,
                   conf = 0.999,
                   max_iters = 50000,
                   laf_consistensy_coef = -1.0,
                   error_type = "sampson",
                   symmetric_error_check = True):
    pts1 = convert_and_check(pts1_)
    pts2 = convert_and_check(pts2_)
    n, dim = pts1.shape
    n2, dim2 = pts2.shape
    assert (n == n2) and (dim == dim2)
    if dim == 2 and laf_consistensy_coef > 0:
        warnings.warn('You set laf_consistensy_coef, but provided only (x,y) keypoints. Skipping LAF check')
        laf_consistensy_coef = 0
    try:
        error_type_int = error_type_dict_homography[error_type.lower()]
    except:
        raise ValueError("Error type should be on of {}. Got {} instead".format(list(error_type_dict_homography.keys()),
                                                                            error_type))
    laf_consistensy_coef = max(0, laf_consistensy_coef)
    H, mask = pydegensac.findHomography_(pts1,
                             pts2,
                             px_th,
                             conf,
                             max_iters,
                             error_type_int,
                             symmetric_error_check,
                             laf_consistensy_coef);
    if np.abs(H).sum() == 0:
        # If we haven`t found any good model, output zeros
        mask = [False]*len(mask)
        return H, mask
    H_out = np.linalg.inv(H.T) 
    return H_out, mask

def findFundamentalMatrix(pts1_,
                   pts2_,
                   px_th = 0.5,
                   conf = 0.9999,
                   max_iters = 100000,
                   laf_consistensy_coef = -1.0,
                   error_type = "sampson",
                   symmetric_error_check = True,
                   enable_degeneracy_check = True):
    pts1 = convert_and_check(pts1_)
    pts2 = convert_and_check(pts2_)
    n, dim = pts1.shape
    n2, dim2 = pts2.shape
    assert (n == n2) and (dim == dim2)
    if dim == 2 and laf_consistensy_coef > 0:
        warnings.warn('You set laf_consistensy_coef, but provided only (x,y) keypoints. Skipping LAF check')
        laf_consistensy_coef = 0
    try:
        error_type_int = error_type_dict_fundamental[error_type.lower()]
    except:
        raise ValueError("Error type should be on of {}. Got {} instead".format(list(error_type_dict_fundamental.keys()),
                                                                            error_type))
    laf_consistensy_coef = max(0, laf_consistensy_coef)
    F, mask = pydegensac.findFundamentalMatrix_(pts1,
                         pts2,
                         px_th,
                         conf,
                         max_iters,
                         error_type_int,
                         symmetric_error_check,
                         laf_consistensy_coef,
                         enable_degeneracy_check);
    if np.abs(F).sum() == 0:
        # If we haven`t found any good model, output zeros
        mask = [False]*n
    return F, mask

