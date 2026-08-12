#include <stdexcept>
#include <string>
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

namespace {

// Input arrays converted to the layout the C core expects, plus the C-side
// scratch/output buffers shared by both estimators. The destructor owns all
// cleanup, so every exit path (including exceptions) releases the memory.
struct ConvertedInput {
    size_t num_tents = 0;
    double *u2 = nullptr;
    double *u2_p1 = nullptr;
    double *u2_p2 = nullptr;
    unsigned char *inl = nullptr;
    int *data_out = nullptr;

    ConvertedInput() = default;
    ConvertedInput(const ConvertedInput&) = delete;
    ConvertedInput& operator=(const ConvertedInput&) = delete;
    ~ConvertedInput() {
        free(data_out);
        delete [] u2;
        delete [] u2_p1;
        delete [] u2_p2;
        delete [] inl;
    }
};

void validate_input(const py::buffer_info &buf1, const py::buffer_info &buf1a,
                    double laf_coef, size_t min_pts) {
    const std::string n_ge = "n>=" + std::to_string(min_pts);

    if ((buf1.ndim != 2) || (buf1a.ndim != 2)) {
        throw std::invalid_argument( "x1y1 and x2y2 must be 2-D arrays with dims [n,2] or [n,6]" );
    }

    size_t NUM_TENTS = buf1.shape[0];
    size_t DIM = buf1.shape[1];

    if ((DIM != 2) && (DIM != 6)) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], [n,6], " + n_ge );
    }
    if (NUM_TENTS < min_pts) {
        throw std::invalid_argument( "x1y1 should be an array with dims [n,2], " + n_ge );
    }
    size_t NUM_TENTSa = buf1a.shape[0];
    size_t DIMa = buf1a.shape[1];

    if ((DIMa != 2) && (DIMa != 6)) {
        throw std::invalid_argument( "x2y2 should be an array with dims [n,2] or [n, 6], " + n_ge );
    }
    if (NUM_TENTSa != NUM_TENTS) {
        throw std::invalid_argument( "x1y1 and x2y2 should be the same size");
    }
    if (DIM != DIMa) {
        throw std::invalid_argument( "x1y1 and x2y2 must have the same number of columns");
    }
    if ((laf_coef > 0) && (DIM == 2)) {
        throw std::invalid_argument( "laf_coef > 0 requires [n,6] input (LAF data)");
    }
}

// Builds the homogeneous [x1 y1 1 x2 y2 1] array (and, with LAFs, the two
// affine-frame-shifted variants) from the Nx2/Nx6 numpy inputs.
void convert_input(const py::buffer_info &buf1, const py::buffer_info &buf1a,
                   double laf_coef, ConvertedInput &in) {
    size_t NUM_TENTS = buf1.shape[0];
    size_t DIM = buf1.shape[1];
    double *ptr1 = (double *) buf1.ptr; // pointer to x1y1 data
    double *ptr1a = (double *) buf1a.ptr; // pointer to x2y2 data

    in.num_tents = NUM_TENTS;

    double *u2Ptr = new double[NUM_TENTS*6];
    in.u2 = u2Ptr;

    // Allocate space only if needed
    int do_laf_check = laf_coef > 0;
    double *u2Ptr_p1 = do_laf_check ? new double[NUM_TENTS*6] : nullptr;
    in.u2_p1 = u2Ptr_p1;
    double *u2Ptr_p2 = do_laf_check ? new double[NUM_TENTS*6] : nullptr;
    in.u2_p2 = u2Ptr_p2;

    typedef unsigned char uchar;
    in.inl = new uchar[NUM_TENTS];

    if (do_laf_check) {
        for (size_t i=0; i < NUM_TENTS; i++) {

            //x1,y1,1
            *u2Ptr =  ptr1[DIM*i];
            u2Ptr++;
            *u2Ptr =  ptr1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x2,y2,1
            *u2Ptr =  ptr1a[DIM*i];
            u2Ptr++;
            *u2Ptr =  ptr1a[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            //x1 + a12,y1 + a22,1
            *u2Ptr_p1 = ptr1[DIM*i] + ptr1[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = ptr1[DIM*i+1] + ptr1[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;

            //x2 + a12,y2 + a22,1
            *u2Ptr_p1 = ptr1a[DIM*i] + ptr1a[DIM*i+3];
            u2Ptr_p1++;
            *u2Ptr_p1 = ptr1a[DIM*i+1] + ptr1a[DIM*i+5];
            u2Ptr_p1++;
            *u2Ptr_p1 =  1.;
            u2Ptr_p1++;


            //x1 + a11,y1 + a21,1
            *u2Ptr_p2 = ptr1[DIM*i] + ptr1[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = ptr1[DIM*i+1] + ptr1[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

            //x2 + a11,y2 + a21,1
            *u2Ptr_p2 = ptr1a[DIM*i] + ptr1a[DIM*i+2];
            u2Ptr_p2++;
            *u2Ptr_p2 = ptr1a[DIM*i+1] + ptr1a[DIM*i+4];
            u2Ptr_p2++;
            *u2Ptr_p2 =  1.;
            u2Ptr_p2++;

        }
    } else {
        for (size_t i=0; i < NUM_TENTS; i++) {

            *u2Ptr =  ptr1[DIM*i];
            u2Ptr++;

            *u2Ptr =  ptr1[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;

            *u2Ptr =  ptr1a[DIM*i];
            u2Ptr++;

            *u2Ptr =  ptr1a[DIM*i+1];
            u2Ptr++;
            *u2Ptr =  1.;
            u2Ptr++;
        };
    }

    in.data_out = (int *) malloc(NUM_TENTS * 18 * sizeof(int));
    if (in.data_out == nullptr) {
        throw std::bad_alloc();
    }
}

py::tuple pack_output(const double *model, const unsigned char *inl, size_t num_tents) {
    //Model
    py::array_t<double> model_out = py::array_t<double>({3,3});
    py::buffer_info buf_model_out = model_out.request();
    double *ptr_model_out = (double *)buf_model_out.ptr;

    for (size_t i=0; i<9; i++)
        ptr_model_out[i]=model[i];

    //Inliers
    py::array_t<bool> inliers_out = py::array_t<bool>(num_tents);
    py::buffer_info buf_inliers = inliers_out.request();
    bool *ptr_inliers= (bool *)buf_inliers.ptr;
    for (size_t i = 0; i < num_tents; i++)
        ptr_inliers[i] = (bool) inl[i];

    return py::make_tuple(model_out, inliers_out);
}

} // namespace

py::tuple findHomography_(py::array_t<double, py::array::c_style | py::array::forcecast>  x1y1_,
                          py::array_t<double, py::array::c_style | py::array::forcecast>   x2y2_,
                          double px_th,
                          double conf,
                          int max_iters,
                          int error_type,
                          bool sym_check_enable,
                          double laf_coef,
                          int seed) {
    // Get the data
    py::buffer_info buf1 = x1y1_.request();
    py::buffer_info buf1a = x2y2_.request();

    validate_input(buf1, buf1a, laf_coef, 4);

    int oriented_constr = 1;
    HDsPtr HDS1 = nullptr;
    HDsiPtr HDSi1 = nullptr;
    HDsidxPtr HDSidx1 = nullptr;

    double error_threshold = 0.0;
    double SymCheck_th = 0.0;
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
    default: {
        throw std::invalid_argument("Unsupported error_type");
    }
    }

    double H[3*3] = {0};

    ConvertedInput in;
    convert_input(buf1, buf1a, laf_coef, in);

    // Run the RANSAC
    exp_ransacHcustomLAF(in.u2,
                         in.u2_p1,
                         in.u2_p2,
                         in.num_tents,
                         error_threshold,
                         laf_coef,
                         conf,
                         max_iters,
                         H,
                         in.inl,
                         4,
                         in.data_out,
                         oriented_constr,
                         0,
                         NULL,
                         HDS1,HDSi1,HDSidx1,
                         SymCheck_th,
                         seed);

    return pack_output(H, in.inl, in.num_tents);
}

py::tuple findFundamentalMatrix_(py::array_t<double, py::array::c_style | py::array::forcecast>  x1y1_,
                                 py::array_t<double, py::array::c_style | py::array::forcecast>  x2y2_,
                                 double px_th,
                                 double conf,
                                 int max_iters,
                                 int error_type,
                                 bool sym_check_enable,
                                 double laf_coef,
                                 bool enable_degeneracy_check,
                                 int seed) {
    // Get the data
    py::buffer_info buf1 = x1y1_.request();
    py::buffer_info buf1a = x2y2_.request();

    validate_input(buf1, buf1a, laf_coef, 8);

    FDsPtr FDS1 = nullptr;
    exFDsPtr EXFDS1 = nullptr;
    FDsidxPtr FDSidx1 = nullptr;
    // SoA form of FDS1, used for the whole-array error evaluations in the
    // RANSAC loop; identical results, contiguous loads (see Ftools.c).
    FDsPtr FDS1soa = nullptr;

    double error_threshold = 0.0;
    double SymCheck_th = 0.0;
    const double SYM_CHECK_COEF = 3.0*sym_check_enable;
    switch (error_type)   {
    case SAMPSON_F:   {
        FDS1 = &FDs;
        EXFDS1 = &exFDs;
        FDSidx1 = &FDsidx;
        FDS1soa = &FDs_soa;

        error_threshold = px_th*px_th;
        SymCheck_th = px_th*px_th * SYM_CHECK_COEF;
        break;
    }

    case SYMM_EPI_F:   {
        FDS1 = &FDsSym;
        EXFDS1 = &exFDsSym;
        FDSidx1 = &FDsSymidx;
        FDS1soa = &FDsSym_soa;
        error_threshold = px_th*px_th;
        SymCheck_th = px_th*px_th * SYM_CHECK_COEF;
        break;
    }
    default: {
        throw std::invalid_argument("Unsupported error_type");
    }
    }

    double F[3*3] = {0};

    ConvertedInput in;
    convert_input(buf1, buf1a, laf_coef, in);

    // Run the RANSAC
    exp_ransacFcustomLAF(in.u2,
                         in.u2_p1,
                         in.u2_p2,
                         in.num_tents,
                         error_threshold,
                         laf_coef,
                         conf,
                         max_iters,
                         F,
                         in.inl,
                         in.data_out,
                         1, 0,
                         NULL,
                         EXFDS1,FDS1,FDSidx1,
                         SymCheck_th,
                         (int)enable_degeneracy_check,
                         seed,
                         FDS1soa);

    return pack_output(F, in.inl, in.num_tents);
}


PYBIND11_MODULE(pydegensac, m) {
    m.doc() = R"doc(
                 Python module
                 -----------------------
                 .. currentmodule:: pydegensac
                 .. autosummary::
                 :toctree: _generate

                 findHomography_,
                 findFundamentalMatrix_

                 )doc";


    m.def("findHomography_", &findHomography_, R"doc(some doc)doc",
            py::arg("x1y1"),
            py::arg("x2y2"),
            py::arg("px_th") = 1.0,
            py::arg("conf") = 0.999,
            py::arg("max_iters") = 10000,
            py::arg("error_type") = 0,
            py::arg("sym_check_enable") = 1,
            py::arg("laf_coef") = 0,
            py::arg("seed") = -1);

    m.def("findFundamentalMatrix_", &findFundamentalMatrix_, R"doc(some doc)doc",
          py::arg("x1y1"),
          py::arg("x2y2"),
          py::arg("px_th") = 0.5,
          py::arg("conf") = 0.9999,
          py::arg("max_iters") = 200000,
          py::arg("error_type") = 0,
          py::arg("sym_check_enable") = 1,
          py::arg("laf_coef") = 0,
          py::arg("enable_degeneracy_check") = 1,
          py::arg("seed") = -1);
}
