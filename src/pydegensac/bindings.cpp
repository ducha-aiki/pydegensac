#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "degensac/exp_ranF.h"
#include "degensac/exp_ranH.h"

namespace py = pybind11;

enum RANSAC_error_t_h {SAMPSON = 0,
                     SYMM_SQ_MAX = 1,
                     SYMM_MAX = 2,
                     SYMM_SQ_SUM = 3,
                     SYMM_SUM = 4};

enum RANSAC_error_t_f {SAMPSON_F = 0,
    SYMM_EPI_F = 1};

py::tuple findHomography_(py::array_t<double>  x1y1_,
                          py::array_t<double>   x2y2_,
                          double px_th,
                          double conf,
                          int max_iters,
                          int error_type,
                          bool sym_check_enable,
                          double laf_coef) {
    // Get the data
    py::buffer_info buf1 = x1y1_.request();
    size_t NUM_TENTS = buf1.shape[0];
    size_t DIM = buf1.shape[1];

    if ((DIM != 2) && (DIM != 6)) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], [n,6], n>=4" );
    }
    if (NUM_TENTS < 4) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], n>=4");
    }
    py::buffer_info buf1a = x2y2_.request();
    size_t NUM_TENTSa = buf1a.shape[0];
    size_t DIMa = buf1a.shape[1];

    if ((DIMa != 2) && (DIMa != 6)) {
        throw std::invalid_argument( "x2y2 should be an array with dims [n,2] or [n, 6], n>=4" );
    }
    if (NUM_TENTSa != NUM_TENTS) {
        throw std::invalid_argument( "x1y1 and x2y2 should be the same size");
    }

    double *ptr1 = (double *) buf1.ptr;
    std::vector<double> x1y1;
    x1y1.assign(ptr1, ptr1 + buf1.size);

    double *ptr1a = (double *) buf1a.ptr;
    std::vector<double> x2y2;
    x2y2.assign(ptr1a, ptr1a + buf1a.size);

    // Convert the data

    int oriented_constr = 1;
    HDsPtr HDS1;
    HDsiPtr HDSi1;
    HDsidxPtr HDSidx1;

    double error_threshold, SymCheck_th;
    const double SYM_CHECK_COEF = 3.0*sym_check_enable;
    switch (error_type)   {
    case SAMPSON:   {
        HDS1 = &HDs;
        HDSi1 = &HDsi;
        HDSidx1 = &HDsidx;
        error_threshold = px_th*px_th;
        SymCheck_th = px_th * SYM_CHECK_COEF;
        break;
    }
    case SYMM_SQ_MAX:   {
        HDS1 = &HDsSymMaxSq;
        HDSi1 = &HDsiSymMaxSq;
        HDSidx1 = &HDsSymMaxSqidx;
        error_threshold = px_th*px_th;
        SymCheck_th = 0;
        break;
    }
    case SYMM_MAX:   {
        HDS1 = &HDsSymMax;
        HDSi1 = &HDsiSymMax;
        HDSidx1 = &HDsSymMaxidx;
        error_threshold = px_th;
        SymCheck_th = 0;
        break;
    }
    case SYMM_SQ_SUM:   {
        HDS1 = &HDsSymSumSq;
        HDSi1 = &HDsiSymSumSq;
        HDSidx1 = &HDsSymSumSqidx;
        error_threshold = px_th*px_th;
        SymCheck_th = px_th * SYM_CHECK_COEF;
        break;
    }
    case SYMM_SUM:   {
        HDS1 = &HDsSymSum;
        HDSi1 = &HDsiSymSum;
        HDSidx1 = &HDsSymSumidx;
        error_threshold = px_th;
        SymCheck_th = px_th * SYM_CHECK_COEF;
        break;
    }
    }


    double H[3*3];

    double *u2Ptr = new double[NUM_TENTS*6], *u2;
    u2=u2Ptr;

    // Allocate space only if needed
    int do_laf_check = laf_coef > 0;
    double *u2Ptr_p1 = new double[do_laf_check*NUM_TENTS*6], *u2_p1;
    u2_p1=u2Ptr_p1;
    double *u2Ptr_p2 = new double[do_laf_check*NUM_TENTS*6], *u2_p2;
    u2_p2=u2Ptr_p2;

    typedef unsigned char uchar;
    unsigned char *inl = new uchar[NUM_TENTS];


    if (do_laf_check) {
        for (size_t i=0; i < NUM_TENTS; i++) {

            //x1,y1,1
            *u2Ptr =  x1y1[DIM*i];
            u2Ptr++;
            *u2Ptr =  x1y1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x2,y2,1
            *u2Ptr =  x2y2[DIM*i];
            u2Ptr++;
            *u2Ptr =  x2y2[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x1 + a12,y1 + a22,1
            *u2Ptr_p1 = x1y1[DIM*i] + x1y1[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = x1y1[DIM*i+1] + x1y1[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;

            //x2 + a12,y2 + a22,1
            *u2Ptr_p1 = x2y2[DIM*i] + x2y2[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = x2y2[DIM*i+1] + x2y2[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;


            //x1 + a11,y1 + a21,1
            *u2Ptr_p2 = x1y1[DIM*i] + x1y1[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = x1y1[DIM*i+1] + x1y1[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

            //x2 + a11,y2 + a21,1
            *u2Ptr_p2 = x2y2[DIM*i] + x2y2[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = x2y2[DIM*i+1] + x2y2[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

        }
    } else {
        for (size_t i=0; i < NUM_TENTS; i++) {

            *u2Ptr =  x1y1[DIM*i];
            u2Ptr++;

            *u2Ptr =  x1y1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            *u2Ptr =  x2y2[DIM*i];
            u2Ptr++;

            *u2Ptr =  x2y2[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;
        };
    }


    int* data_out = (int *) malloc(NUM_TENTS * 18 * sizeof(int));
    double *resids;


    // Run the RANSAC
    exp_ransacHcustomLAF(u2,
                         u2_p1,
                         u2_p2,
                         NUM_TENTS,
                         error_threshold,
                         laf_coef,
                         conf,
                         max_iters,
                         H,
                         inl,
                         4,
                         data_out,
                         oriented_constr,
                         0,
                         &resids,
                         HDS1,HDSi1,HDSidx1,
                         SymCheck_th);



    //Model
    py::array_t<double> H_out = py::array_t<double>({3,3});
    py::buffer_info buf_H_out = H_out.request();
    double *ptr_H_out = (double *)buf_H_out.ptr;

    for (size_t i=0; i<9; i++)
        ptr_H_out[i]=H[i];

    //Inliers
    py::array_t<bool> inliers_out = py::array_t<bool>(NUM_TENTS);
    py::buffer_info buf_inliers = inliers_out.request();
    bool *ptr_inliers= (bool *)buf_inliers.ptr;
    for (size_t i = 0; i < NUM_TENTS; i++)
        ptr_inliers[i] = (bool) inl[i];


    free(resids);
    free(data_out);
    delete [] u2;
    delete [] u2_p1;
    delete [] u2_p2;
    delete [] inl;


    return py::make_tuple(H_out, inliers_out);
}

py::tuple findFundamentalMatrix_(py::array_t<double>  x1y1_,
                                 py::array_t<double>   x2y2_,
                                 double px_th,
                                 double conf,
                                 int max_iters,
                                 int error_type,
                                 bool sym_check_enable,
                                 double laf_coef,
                                 bool enable_degeneracy_check) {
    // Get the data
    py::buffer_info buf1 = x1y1_.request();
    size_t NUM_TENTS = buf1.shape[0];
    size_t DIM = buf1.shape[1];

    if ((DIM != 2) && (DIM != 6)) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], [n,6], n>=8" );
    }
    if (NUM_TENTS < 8) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], n>=8");
    }
    py::buffer_info buf1a = x2y2_.request();
    size_t NUM_TENTSa = buf1a.shape[0];
    size_t DIMa = buf1a.shape[1];

    if ((DIMa != 2) && (DIMa != 6)) {
        throw std::invalid_argument( "x2y2 should be an array with dims [n,2] or [n, 6], n>=8" );
    }
    if (NUM_TENTSa != NUM_TENTS) {
        throw std::invalid_argument( "x1y1 and x2y2 should be the same size");
    }

    double *ptr1 = (double *) buf1.ptr;
    std::vector<double> x1y1;
    x1y1.assign(ptr1, ptr1 + buf1.size);

    double *ptr1a = (double *) buf1a.ptr;
    std::vector<double> x2y2;
    x2y2.assign(ptr1a, ptr1a + buf1a.size);

    // Convert the data
    FDsPtr FDS1;
    exFDsPtr EXFDS1;
    FDsidxPtr FDSidx1;

    double error_threshold, SymCheck_th;
    const double SYM_CHECK_COEF = 3.0*sym_check_enable;
    switch (error_type)   {
    case SAMPSON_F:   {
        FDS1 = &FDs;
        EXFDS1 = &exFDs;
        FDSidx1 = &FDsidx;

        error_threshold = px_th*px_th;
        SymCheck_th = px_th*px_th * SYM_CHECK_COEF;
        break;
    }

    case SYMM_EPI_F:   {
        FDS1 = &FDsSym;
        EXFDS1 = &exFDsSym;
        FDSidx1 = &FDsSymidx;
        error_threshold = px_th*px_th;
        SymCheck_th = px_th*px_th * SYM_CHECK_COEF;
        break;
    }
    }


    double F[3*3];

    double *u2Ptr = new double[NUM_TENTS*6], *u2;
    u2=u2Ptr;

    // Allocate space only if needed
    int do_laf_check = laf_coef > 0;
    double *u2Ptr_p1 = new double[do_laf_check*NUM_TENTS*6], *u2_p1;
    u2_p1=u2Ptr_p1;
    double *u2Ptr_p2 = new double[do_laf_check*NUM_TENTS*6], *u2_p2;
    u2_p2=u2Ptr_p2;

    typedef unsigned char uchar;
    unsigned char *inl = new uchar[NUM_TENTS];


    if (do_laf_check) {
        for (size_t i=0; i < NUM_TENTS; i++) {

            //x1,y1,1
            *u2Ptr =  x1y1[DIM*i];
            u2Ptr++;
            *u2Ptr =  x1y1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x2,y2,1
            *u2Ptr =  x2y2[DIM*i];
            u2Ptr++;
            *u2Ptr =  x2y2[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x1 + a12,y1 + a22,1
            *u2Ptr_p1 = x1y1[DIM*i] + x1y1[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = x1y1[DIM*i+1] + x1y1[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;

            //x2 + a12,y2 + a22,1
            *u2Ptr_p1 = x2y2[DIM*i] + x2y2[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = x2y2[DIM*i+1] + x2y2[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;


            //x1 + a11,y1 + a21,1
            *u2Ptr_p2 = x1y1[DIM*i] + x1y1[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = x1y1[DIM*i+1] + x1y1[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

            //x2 + a11,y2 + a21,1
            *u2Ptr_p2 = x2y2[DIM*i] + x2y2[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = x2y2[DIM*i+1] + x2y2[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

        }
    } else {
        for (size_t i=0; i < NUM_TENTS; i++) {

            *u2Ptr =  x1y1[DIM*i];
            u2Ptr++;

            *u2Ptr =  x1y1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            *u2Ptr =  x2y2[DIM*i];
            u2Ptr++;

            *u2Ptr =  x2y2[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;
        };
    }


    int* data_out = (int *) malloc(NUM_TENTS * 18 * sizeof(int));
    double *resids;
    int I_H = 0;
    int *Ihptr = &I_H;
    double HinF [3*3];

    // Run the RANSAC
    exp_ransacFcustomLAF(u2,
                         u2_p1,
                         u2_p2,
                         NUM_TENTS,
                         error_threshold,
                         laf_coef,
                         conf,
                         max_iters,
                         F,
                         inl,
                         data_out,
                         1, 0,
                         &resids,
                         HinF,Ihptr,
                         EXFDS1,FDS1,FDSidx1,
                         SymCheck_th,
                         (int)enable_degeneracy_check);


    // Convert and store output



    //Model
    py::array_t<double> F_out = py::array_t<double>({3,3});
    py::buffer_info buf_F_out = F_out.request();
    double *ptr_F_out = (double *)buf_F_out.ptr;

    for (size_t i=0; i<9; i++)
        ptr_F_out[i]=F[i];

    //Inliers
    py::array_t<bool> inliers_out = py::array_t<bool>(NUM_TENTS);
    py::buffer_info buf_inliers = inliers_out.request();
    bool *ptr_inliers= (bool *)buf_inliers.ptr;
    for (size_t i = 0; i < NUM_TENTS; i++)
        ptr_inliers[i] = (bool) inl[i];


    free(resids);
    free(data_out);
    delete [] u2;
    delete [] u2_p1;
    delete [] u2_p2;
    delete [] inl;


    return py::make_tuple(F_out, inliers_out);
}


PYBIND11_PLUGIN(pydegensac) {
    py::module m("pydegensac", R"doc(
                 Python module
                 -----------------------
                 .. currentmodule:: pydegensac
                 .. autosummary::
                 :toctree: _generate

                 findHomography_,
                 findFundamentalMatrix_

                 )doc");


    m.def("findHomography_", &findHomography_, R"doc(some doc)doc",
            py::arg("x1y1"),
            py::arg("x2y2"),
            py::arg("px_th") = 1.0,
            py::arg("conf") = 0.999,
            py::arg("max_iters") = 10000,
            py::arg("error_type") = 0,
            py::arg("sym_check_enable") = 1,
            py::arg("laf_coef") = 0);

    m.def("findFundamentalMatrix_", &findFundamentalMatrix_, R"doc(some doc)doc",
          py::arg("x1y1"),
          py::arg("x2y2"),
          py::arg("px_th") = 0.5,
          py::arg("conf") = 0.9999,
          py::arg("max_iters") = 200000,
          py::arg("error_type") = 0,
          py::arg("sym_check_enable") = 1,
          py::arg("laf_coef") = 0,
          py::arg("enable_degeneracy_check") = 1);

    return m.ptr();
}
