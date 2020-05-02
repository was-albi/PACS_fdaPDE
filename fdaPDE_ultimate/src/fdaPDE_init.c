#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP eval_FEM_fd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP points_projection(SEXP , SEXP);
extern SEXP tree_mesh_construction(SEXP , SEXP , SEXP , SEXP );
extern SEXP get_FEM_mass_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_space_varying_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_stiff_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_integration_points(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_triangulate_native(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_Laplace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE_space_varying(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Smooth_FPCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_Laplace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_PDE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gam_PDE_space_varying( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_regression_Estimates(SEXP , SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"eval_FEM_fd",                      (DL_FUNC) &eval_FEM_fd,                       8},
    {"points_projection",                (DL_FUNC) &points_projection,                 2},
    {"tree_mesh_construction",           (DL_FUNC) &tree_mesh_construction,            4},
    {"get_FEM_mass_matrix",              (DL_FUNC) &get_FEM_mass_matrix,               4},
    {"get_FEM_PDE_matrix",               (DL_FUNC) &get_FEM_PDE_matrix,               17},
    {"get_FEM_PDE_space_varying_matrix", (DL_FUNC) &get_FEM_PDE_space_varying_matrix, 18},
    {"get_FEM_stiff_matrix",             (DL_FUNC) &get_FEM_stiff_matrix,              4},
    {"get_integration_points",           (DL_FUNC) &get_integration_points,            4},
    {"R_triangulate_native",             (DL_FUNC) &R_triangulate_native,              8},
    {"regression_Laplace",               (DL_FUNC) &regression_Laplace,               14},
    {"regression_PDE",                   (DL_FUNC) &regression_PDE,                   17},
    {"regression_PDE_space_varying",     (DL_FUNC) &regression_PDE_space_varying,     18},
    {"Smooth_FPCA",                      (DL_FUNC) &Smooth_FPCA,                      13},
    {"gam_Laplace",                      (DL_FUNC) &gam_Laplace,                      20},
    {"gam_PDE",                          (DL_FUNC) &gam_PDE,                          23},
    {"gam_PDE_space_varying",            (DL_FUNC) &gam_PDE_space_varying,            24},
    {"get_regression_Estimates",         (DL_FUNC) &get_regression_Estimates,         10},
    {NULL, NULL, 0}
};

void R_init_fdaPDE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
