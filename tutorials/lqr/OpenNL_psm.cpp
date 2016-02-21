#include "OpenNL_psm.h"

/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */


/*
 *  This file is a PSM (pluggable software module)
 *   generated from the distribution of Geogram.
 *
 *  See Geogram documentation on:
 *   http://alice.loria.fr/software/geogram/doc/html/index.html
 *
 *  See documentation of the functions bundled in this PSM on:
 *   http://alice.loria.fr/software/geogram/doc/html/nl_8h.html
 */



/******* extracted from nl_private.h *******/

#ifndef __NL_PRIVATE__
#define __NL_PRIVATE__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>



#define nl_arg_used(x) (void)x

#define nl_cast(T,x)   (*(T*)&(x))



void nl_assertion_failed(const char* cond, const char* file, int line) ;

void nl_range_assertion_failed(
    double x, double min_val, double max_val, const char* file, int line
) ;

void nl_should_not_have_reached(const char* file, int line) ;

#define nl_assert(x) {                                          \
    if(!(x)) {                                                  \
        nl_assertion_failed(#x,__FILE__, __LINE__) ;            \
    }                                                           \
} 

#define nl_range_assert(x,min_val,max_val) {                    \
    if(((x) < (min_val)) || ((x) > (max_val))) {                \
        nl_range_assertion_failed(x, min_val, max_val,          \
            __FILE__, __LINE__                                  \
        ) ;                                                     \
    }                                                           \
}

#define nl_assert_not_reached {                                 \
    nl_should_not_have_reached(__FILE__, __LINE__) ;            \
}

#ifdef NL_DEBUG
    #define nl_debug_assert(x) nl_assert(x)
    #define nl_debug_range_assert(x,min_val,max_val)            \
                               nl_range_assert(x,min_val,max_val)
#else
    #define nl_debug_assert(x) 
    #define nl_debug_range_assert(x,min_val,max_val) 
#endif

#ifdef NL_PARANOID
    #define nl_parano_assert(x) nl_assert(x)
    #define nl_parano_range_assert(x,min_val,max_val)           \
                               nl_range_assert(x,min_val,max_val)
#else
    #define nl_parano_assert(x) 
    #define nl_parano_range_assert(x,min_val,max_val) 
#endif


void nlError(const char* function, const char* message) ;

void nlWarning(const char* function, const char* message) ;


NLdouble nlCurrentTime()  ;


/* classic macros */

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y)) 
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#endif



#define NL_NEW(T)                (T*)(calloc(1, sizeof(T))) 

#define NL_NEW_ARRAY(T,NB)       (T*)(calloc((size_t)(NB),sizeof(T)))

#define NL_RENEW_ARRAY(T,x,NB)   (T*)(realloc(x,(size_t)(NB)*sizeof(T))) 

#define NL_DELETE(x)             free(x); x = NULL 

#define NL_DELETE_ARRAY(x)       free(x); x = NULL

#define NL_CLEAR(T, x)           memset(x, 0, sizeof(T)) 

#define NL_CLEAR_ARRAY(T,x,NB)   memset(x, 0, (size_t)(NB)*sizeof(T)) 



#define NL_UINT_MAX 0xffffffff

#define NL_USHORT_MAX 0xffff


#endif

/******* extracted from nl_matrix.h *******/


#ifndef __NL_MATRIX__
#define __NL_MATRIX__


#ifdef __cplusplus
extern "C" {
#endif


/* Dynamic arrays for sparse row/columns */

typedef struct  {
    NLuint index ;

    NLdouble value ; 
} NLCoeff ;

typedef struct {
    NLuint size ;
    
    NLuint capacity ;

    NLCoeff* coeff ;  
} NLRowColumn ;

void nlRowColumnConstruct(NLRowColumn* c) ;

void nlRowColumnDestroy(NLRowColumn* c) ;

void nlRowColumnGrow(NLRowColumn* c) ;

void nlRowColumnAdd(NLRowColumn* c, NLuint index, NLdouble value) ;

void nlRowColumnAppend(NLRowColumn* c, NLuint index, NLdouble value) ;

void nlRowColumnZero(NLRowColumn* c) ;

void nlRowColumnClear(NLRowColumn* c) ;

void nlRowColumnSort(NLRowColumn* c) ;


/* Compressed Row Storage */

typedef struct {
    NLuint m;
    
    NLuint n;

    NLdouble* val;    

    NLuint* rowptr;

    NLuint* colind;

    NLuint nslices;

    NLuint* sliceptr; 
} NLCRSMatrix;

void nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
);

void nlCRSMatrixDestroy(NLCRSMatrix* M);

NLboolean nlCRSMatrixLoad(NLCRSMatrix* M, const char* filename);

NLboolean nlCRSMatrixSave(NLCRSMatrix* M, const char* filename);


/* SparseMatrix data structure */

#define NL_MATRIX_STORE_ROWS          1

#define NL_MATRIX_STORE_COLUMNS       2

#define NL_MATRIX_STORE_SYMMETRIC     4

#define NL_MATRIX_STORE_COMPRESSED    8

#define NL_MATRIX_STORE_DIAG_INV      16
    
typedef struct {
    NLuint m ;
    
    NLuint n ;

    NLuint diag_size ;

    NLenum storage ;

    NLRowColumn* row ;

    NLRowColumn* column ;

    NLdouble*    diag ;

    NLdouble*    diag_inv ;

    NLCRSMatrix* compressed ;  
} NLSparseMatrix ;


void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) ;

void nlSparseMatrixDestroy(NLSparseMatrix* M) ;

void nlSparseMatrixAdd(
    NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value
) ;

void nlSparseMatrixZero( NLSparseMatrix* M) ;

void nlSparseMatrixClear( NLSparseMatrix* M) ;

NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) ;

void nlSparseMatrixSort( NLSparseMatrix* M) ;

void nlSparseMatrixComputeDiagInv( NLSparseMatrix* M);

void nlSparseMatrixCompress( NLSparseMatrix* M);
    

/* SparseMatrix x Vector routine */

void nlSparseMatrixMult(NLSparseMatrix* A, const NLdouble* x, NLdouble* y) ;

#ifdef __cplusplus
}
#endif

#endif

/******* extracted from nl_context.h *******/

#ifndef __NL_CONTEXT__
#define __NL_CONTEXT__




/* NLContext data structure */

typedef void(*NLMatrixFunc)(const double* x, double* y);

typedef NLboolean(*NLSolverFunc)();

typedef void(*NLProgressFunc)(
    NLuint cur_iter, NLuint max_iter, double cur_err, double max_err
);

typedef struct {
    NLdouble  value ;
    NLboolean locked ;
    NLuint    index ;
} NLVariable ;

#define NL_STATE_INITIAL                0
#define NL_STATE_SYSTEM                 1
#define NL_STATE_MATRIX                 2
#define NL_STATE_ROW                    3
#define NL_STATE_MATRIX_CONSTRUCTED     4
#define NL_STATE_SYSTEM_CONSTRUCTED     5
#define NL_STATE_SOLVED                 6

typedef struct {
    NLenum           state ;

    NLVariable*      variable ;

    NLuint           n ;

    NLSparseMatrix   M ;

    NLRowColumn      af ;

    NLRowColumn      al ;

    NLRowColumn      xl ;

    NLdouble*        x ;

    NLdouble*        b ;

    NLdouble         right_hand_side ;

    NLboolean        right_hand_side_set ;
    
    NLdouble         row_scaling ;

    NLenum           solver ;

    NLenum           preconditioner ;

    NLboolean        preconditioner_defined;
    
    NLuint           nb_variables ;


    NLboolean        ij_coefficient_called;
    
    NLuint           current_row ;

    NLboolean        least_squares ;

    NLboolean        symmetric ;

    NLuint           max_iterations ;


    NLboolean        max_iterations_defined ;
    
    NLuint           inner_iterations ;

    NLdouble         threshold ;

    NLboolean        threshold_defined ;
    
    NLdouble         omega ;

    NLboolean        normalize_rows ;

    NLboolean        alloc_M ;

    NLboolean        alloc_af ;

    NLboolean        alloc_al ;

    NLboolean        alloc_xl ;

    NLboolean        alloc_variable ;

    NLboolean        alloc_x ;

    NLboolean        alloc_b ;

    NLuint           used_iterations ;

    NLdouble         error ;

    NLdouble         elapsed_time ;

    NLMatrixFunc     matrix_vector_prod ;

    NLMatrixFunc     precond_vector_prod ;

    NLSolverFunc     solver_func ;

    NLProgressFunc   progress_func ;

    NLboolean        verbose;

    NLulong          flops;
} NLContextStruct ;

extern NLContextStruct* nlCurrentContext ;

void nlCheckState(NLenum state) ;

void nlTransition(NLenum from_state, NLenum to_state) ;

void nlMatrixVectorProd_default(const NLdouble* x, NLdouble* y) ;

NLboolean nlDefaultSolver() ;

#endif

/******* extracted from nl_blas.h *******/



#ifndef __NL_BLAS__
#define __NL_BLAS__

#ifndef NL_FORTRAN_WRAP
#define NL_FORTRAN_WRAP(x) x##_
#endif



/* C wrappers for BLAS routines */

void dscal( int n, double a, double *x, int incx ) ;

void dcopy( 
    int n, const double *x, int incx, double *y, int incy 
) ;

void daxpy( 
    int n, double a, const double *x, int incx, double *y,
    int incy 
) ;

double ddot( 
    int n, const double *x, int incx, const double *y, int incy 
) ;

double dnrm2( int n, const double *x, int incx ) ;

typedef enum {
    NoTranspose=0, Transpose=1, ConjugateTranspose=2
} MatrixTranspose ;

typedef enum {
    UpperTriangle=0, LowerTriangle=1
} MatrixTriangle ;


typedef enum {
    UnitTriangular=0, NotUnitTriangular=1
} MatrixUnitTriangular ;


void dtpsv( 
    MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
) ;

void dgemv( 
    MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
) ;

#endif

/******* extracted from nl_iterative_solvers.h *******/


#ifndef __NL_ITERATIVE_SOLVERS__
#define __NL_ITERATIVE_SOLVERS__



NLuint nlSolve_CG() ;

NLuint nlSolve_CG_precond() ;

NLuint nlSolve_BICGSTAB() ;

NLuint nlSolve_BICGSTAB_precond() ;

NLuint nlSolve_GMRES() ;

#endif


/******* extracted from nl_preconditioners.h *******/

#ifndef __NL_PRECONDITIONERS__
#define __NL_PRECONDITIONERS__




/* preconditioners */

void nlPreconditioner_Jacobi(const NLdouble* x, NLdouble* y) ;

void nlPreconditioner_SSOR(const NLdouble* x, NLdouble* y) ;


void nlMultDiagonal(NLdouble* x, NLdouble omega) ;

void nlMultDiagonalInverse(NLdouble* x, NLdouble omega) ;

void nlMultLowerInverse(const NLdouble* x, NLdouble* y, NLdouble omega) ;

void nlMultUpperInverse(const NLdouble* x, NLdouble* y, NLdouble omega) ;

#endif

/******* extracted from nl_superlu.h *******/

#ifndef __NL_SUPERLU__
#define __NL_SUPERLU__



NLboolean nlSolve_SUPERLU();

NLboolean nlInitExtension_SUPERLU();


#endif

/******* extracted from nl_cnc_gpu_cuda.h *******/

#ifndef __NL_CNC_GPU_CUDA__
#define __NL_CNC_GPU_CUDA__

#ifdef __cplusplus
extern "C" {
#endif


    
NLboolean nlSolverIsCNC(NLint solver) ;

NLuint nlSolve_CNC() ;

#ifdef __cplusplus
}
#endif

#endif

/******* extracted from nl_os.c *******/



#if (defined (WIN32) || defined(_WIN64))
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/times.h> 
#endif


/* Assertions */


void nl_assertion_failed(const char* cond, const char* file, int line) {
    fprintf(
        stderr, 
        "OpenNL assertion failed: %s, file:%s, line:%d\n",
        cond,file,line
    ) ;
    abort() ;
}

void nl_range_assertion_failed(
    double x, double min_val, double max_val, const char* file, int line
) {
    fprintf(
        stderr, 
        "OpenNL range assertion failed: %f in [ %f ... %f ], file:%s, line:%d\n",
        x, min_val, max_val, file,line
    ) ;
    abort() ;
}

void nl_should_not_have_reached(const char* file, int line) {
    fprintf(
        stderr, 
        "OpenNL should not have reached this point: file:%s, line:%d\n",
        file,line
    ) ;
    abort() ;
}



/* Timing */

#ifdef WIN32
NLdouble nlCurrentTime() {
    return (NLdouble)GetTickCount() / 1000.0 ;
}
#else
double nlCurrentTime() {
    clock_t user_clock ;
    struct tms user_tms ;
    user_clock = times(&user_tms) ;
    return (NLdouble)user_clock / 100.0 ;
}
#endif



/* Error-reporting functions */

void nlError(const char* function, const char* message) {
    fprintf(stderr, "OpenNL error in %s(): %s\n", function, message) ; 
}

void nlWarning(const char* function, const char* message) {
    fprintf(stderr, "OpenNL warning in %s(): %s\n", function, message) ; 
}





/******* extracted from nl_matrix.c *******/


void nlRowColumnConstruct(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    c->coeff    = NULL ;
}

void nlRowColumnDestroy(NLRowColumn* c) {
    NL_DELETE_ARRAY(c->coeff) ;
#ifdef NL_PARANOID
    NL_CLEAR(NLRowColumn, c) ; 
#endif
}

void nlRowColumnGrow(NLRowColumn* c) {
    if(c->capacity != 0) {
        c->capacity = 2 * c->capacity ;
        c->coeff = NL_RENEW_ARRAY(NLCoeff, c->coeff, c->capacity) ;
    } else {
        c->capacity = 4 ;
        c->coeff = NL_NEW_ARRAY(NLCoeff, c->capacity) ;
    }
}

void nlRowColumnAdd(NLRowColumn* c, NLuint index, NLdouble value) {
    NLuint i ;
    for(i=0; i<c->size; i++) {
        if(c->coeff[i].index == index) {
            c->coeff[i].value += value ;
            return ;
        }
    }
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

/* Does not check whether the index already exists */
void nlRowColumnAppend(NLRowColumn* c, NLuint index, NLdouble value) {
    if(c->size == c->capacity) {
        nlRowColumnGrow(c) ;
    }
    c->coeff[c->size].index = index ;
    c->coeff[c->size].value = value ;
    c->size++ ;
}

void nlRowColumnZero(NLRowColumn* c) {
    c->size = 0 ;
}

void nlRowColumnClear(NLRowColumn* c) {
    c->size     = 0 ;
    c->capacity = 0 ;
    NL_DELETE_ARRAY(c->coeff) ;
}

static int nlCoeffCompare(const void* p1, const void* p2) {
    return (((NLCoeff*)(p2))->index < ((NLCoeff*)(p1))->index) ;
}

void nlRowColumnSort(NLRowColumn* c) {
    qsort(c->coeff, c->size, sizeof(NLCoeff), nlCoeffCompare) ;
}


/* CRSMatrix data structure */

void nlCRSMatrixConstruct(
    NLCRSMatrix* M, NLuint m, NLuint n, NLuint nnz, NLuint nslices
) {
    M->m = m;
    M->n = n;
    M->nslices = nslices;
    M->val = NL_NEW_ARRAY(double, nnz);
    M->rowptr = NL_NEW_ARRAY(NLuint, m+1);
    M->colind = NL_NEW_ARRAY(NLuint, nnz);
    M->sliceptr = NL_NEW_ARRAY(NLuint, nslices+1);
}

void nlCRSMatrixDestroy(NLCRSMatrix* M) {
    NL_DELETE_ARRAY(M->val);
    NL_DELETE_ARRAY(M->rowptr);
    NL_DELETE_ARRAY(M->colind);
    NL_DELETE_ARRAY(M->sliceptr);
    M->m = 0;
    M->n = 0;
    M->nslices = 0;
}


NLboolean nlCRSMatrixSave(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = M->rowptr[M->m];
    FILE* f = fopen(filename, "rb");
    if(f == NULL) {
        nlError("nlCRSMatrixSave", "Could not open file");
        return NL_FALSE;
    }

    fwrite(&M->m, sizeof(NLuint), 1, f);
    fwrite(&M->n, sizeof(NLuint), 1, f);
    fwrite(&nnz, sizeof(NLuint), 1, f);

    fwrite(M->rowptr, sizeof(NLuint), M->m+1, f);
    fwrite(M->colind, sizeof(NLuint), nnz, f);
    fwrite(M->val, sizeof(double), nnz, f);
    
    return NL_TRUE;
}

NLboolean nlCRSMatrixLoad(NLCRSMatrix* M, const char* filename) {
    NLuint nnz = 0;
    FILE* f = fopen(filename, "rb");
    NLboolean truncated = NL_FALSE;
    
    if(f == NULL) {
        nlError("nlCRSMatrixLoad", "Could not open file");
        return NL_FALSE;
    }
    
    truncated = truncated || (
        fread(&M->m, sizeof(NLuint), 1, f) != 1 ||
        fread(&M->n, sizeof(NLuint), 1, f) != 1 ||
        fread(&nnz, sizeof(NLuint), 1, f) != 1
    );

    if(truncated) {
        M->rowptr = NULL;
        M->colind = NULL;
        M->val = NULL;
    } else {
        M->rowptr = NL_NEW_ARRAY(NLuint, M->m+1);
        M->colind = NL_NEW_ARRAY(NLuint, nnz);
        M->val = NL_NEW_ARRAY(double, nnz);
        truncated = truncated || (
            fread(M->rowptr, sizeof(NLuint), M->m+1, f) != M->m+1 ||
            fread(M->colind, sizeof(NLuint), nnz, f) != nnz ||
            fread(M->val, sizeof(double), nnz, f) != nnz
        );
    }

    if(truncated) {
        nlError("nlCRSMatrixSave", "File appears to be truncated");
        NL_DELETE_ARRAY(M->rowptr);
        NL_DELETE_ARRAY(M->colind);
        NL_DELETE_ARRAY(M->val);
        return NL_FALSE;
    } else {
        M->nslices = 1;    
        M->sliceptr = NL_NEW_ARRAY(NLuint, M->nslices+1);
        M->sliceptr[0] = 0;
        M->sliceptr[1] = M->m;
    }

    fclose(f);
    return NL_TRUE;
}

NLuint nlCRSMatrixNNZ(NLCRSMatrix* M) {
    return M->rowptr[M->m];
}

void nlCRSMatrixMultSlice(
    NLCRSMatrix* M, const double* x, double* y, NLuint Ibegin, NLuint Iend
) {
    NLuint i,j;
    for(i=Ibegin; i<Iend; ++i) {
        double sum=0.0;
        for(j=M->rowptr[i]; j<M->rowptr[i+1]; ++j) {
            sum += M->val[j] * x[M->colind[j]];
        }
        y[i] = sum; 
    }
}


void nlCRSMatrixMult(
    NLCRSMatrix* M, const double* x, double* y
) {
    
    int slice;
    int nslices = (int)(M->nslices);
#pragma omp parallel for private(slice) 
    for(slice=0; slice<nslices; ++slice) {
        nlCRSMatrixMultSlice(
            M,x,y,M->sliceptr[slice],M->sliceptr[slice+1]
        );
    }
}


/* SparseMatrix data structure */

void nlSparseMatrixConstruct(
    NLSparseMatrix* M, NLuint m, NLuint n, NLenum storage
) {
    NLuint i ;
    M->m = m ;
    M->n = n ;
    M->storage = storage ;
    if(storage & NL_MATRIX_STORE_ROWS) {
        M->row = NL_NEW_ARRAY(NLRowColumn, m) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->row[i])) ;
        }
    } else {
        M->row = NULL ;
    }

    if(storage & NL_MATRIX_STORE_COLUMNS) {
        M->column = NL_NEW_ARRAY(NLRowColumn, n) ;
        for(i=0; i<n; i++) {
            nlRowColumnConstruct(&(M->column[i])) ;
        }
    } else {
        M->column = NULL ;
    }

    M->diag_size = MIN(m,n) ;
    M->diag = NL_NEW_ARRAY(NLdouble, M->diag_size) ;
}

void nlSparseMatrixDestroyRowColumns(NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnDestroy(&(M->row[i])) ;
        }
        NL_DELETE_ARRAY(M->row) ;
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_ROWS);
    
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnDestroy(&(M->column[i])) ;
        }
        NL_DELETE_ARRAY(M->column) ;
    }
    M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_COLUMNS);    
}

void nlSparseMatrixDestroy(NLSparseMatrix* M) {
    nlSparseMatrixDestroyRowColumns(M);
    NL_DELETE_ARRAY(M->diag) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_DELETE_ARRAY(M->diag_inv) ;
        M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_DIAG_INV);                    
    }    
    if(M->storage & NL_MATRIX_STORE_COMPRESSED) {
        nlCRSMatrixDestroy(M->compressed);
        NL_DELETE(M->compressed);
        M->storage = (NLenum)((int)(M->storage) & ~NL_MATRIX_STORE_COMPRESSED);            
    }
#ifdef NL_PARANOID
    NL_CLEAR(NLSparseMatrix,M) ;
#endif
}

void nlSparseMatrixAdd(NLSparseMatrix* M, NLuint i, NLuint j, NLdouble value) {
    nl_parano_range_assert(i, 0, M->m - 1) ;
    nl_parano_range_assert(j, 0, M->n - 1) ;
    nl_debug_assert(!(M->storage & NL_MATRIX_STORE_COMPRESSED));
    if((M->storage & NL_MATRIX_STORE_SYMMETRIC) && (j > i)) {
        return ;
    }
    if(i == j) {
        M->diag[i] += value ;
    }
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        nlRowColumnAdd(&(M->row[i]), j, value) ;
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        nlRowColumnAdd(&(M->column[j]), i, value) ;
    }
}

void nlSparseMatrixZero( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnZero(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnZero(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_CLEAR_ARRAY(NLdouble, M->diag_inv, M->diag_size) ;        
    }
}

void nlSparseMatrixClear( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i=0; i<M->m; i++) {
            nlRowColumnClear(&(M->row[i])) ;
        }
    }
    if(M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i=0; i<M->n; i++) {
            nlRowColumnClear(&(M->column[i])) ;
        }
    }
    NL_CLEAR_ARRAY(NLdouble, M->diag, M->diag_size) ;
    if(M->storage & NL_MATRIX_STORE_DIAG_INV) {
        NL_CLEAR_ARRAY(NLdouble, M->diag_inv, M->diag_size) ;        
    }
}

/* Returns the number of non-zero coefficients */
NLuint nlSparseMatrixNNZ( NLSparseMatrix* M) {
    NLuint nnz = 0 ;
    NLuint i ;
    if (M->storage & NL_MATRIX_STORE_COMPRESSED) {
        nnz = nlCRSMatrixNNZ(M->compressed);
    } else if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nnz += M->row[i].size ;
        }
    } else if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nnz += M->column[i].size ;
        }
    } else {
        nl_assert_not_reached ;
    }
    return nnz ;
}

void nlSparseMatrixSort( NLSparseMatrix* M) {
    NLuint i ;
    if(M->storage & NL_MATRIX_STORE_ROWS) {
        for(i = 0; i<M->m; i++) {
            nlRowColumnSort(&(M->row[i])) ;                
        }
    } 
    if (M->storage & NL_MATRIX_STORE_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nlRowColumnSort(&(M->column[i])) ;
        }
    } 
}

void nlSparseMatrixComputeDiagInv( NLSparseMatrix* M) {
    NLuint i;
    NLdouble s;
    if(!(M->storage & NL_MATRIX_STORE_DIAG_INV)) {
        M->diag_inv = NL_NEW_ARRAY(double, M->diag_size);
        M->storage |= NL_MATRIX_STORE_DIAG_INV;
    }
    for(i=0; i<M->diag_size; ++i) {
        s = M->diag[i];
        if(s != 0.0) {
            s = 1.0 / s;
        }
        M->diag_inv[i] = s;
    }
}

void nlSparseMatrixCompress( NLSparseMatrix* M) {
    NLuint nnz = nlSparseMatrixNNZ(M);
    NLuint nslices = 8; /* TODO: get number of cores */
    NLuint slice, cur_bound, cur_NNZ, cur_row;
    NLuint i,ij,k;
    NLuint slice_size = nnz / nslices;
    NLCRSMatrix* CRS = NL_NEW(NLCRSMatrix);
    
    nlCRSMatrixConstruct(CRS, M->m, M->n, nnz, nslices);
    
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS);
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC));
    nl_assert(!(M->storage & NL_MATRIX_STORE_COMPRESSED));

    nlSparseMatrixSort(M);

    /* Copy dynamic sparse matrix into compressed row storage */
    k=0;
    for(i=0; i<M->m; ++i) {
        NLRowColumn* Ri = &(M->row[i]) ;
        CRS->rowptr[i] = k;
        for(ij=0; ij<Ri->size; ij++) {
            NLCoeff* c = &(Ri->coeff[ij]) ;
            CRS->val[k] = c->value;
            CRS->colind[k] = c->index;
            ++k;
        }
    }
    CRS->rowptr[M->m] = nnz;

    /* Create "slices" to be used by parallel sparse matrix vector product */
    cur_bound = slice_size;
    cur_NNZ = 0;
    cur_row = 0;
    CRS->sliceptr[0]=0;
    for(slice=1; slice<nslices; ++slice) {
        while(cur_NNZ < cur_bound && cur_row < M->m) {
            ++cur_row;
            cur_NNZ += CRS->rowptr[cur_row+1] - CRS->rowptr[cur_row];
        }
        CRS->sliceptr[slice] = cur_row;
        cur_bound += slice_size;
    }
    CRS->sliceptr[nslices]=M->m; 
    
    M->compressed = CRS;
    M->storage = (M->storage | NL_MATRIX_STORE_COMPRESSED);
    nlSparseMatrixDestroyRowColumns(M);
}




/* SparseMatrix x Vector routines, internal helper routines */

static void nlSparseMatrix_mult_rows_symmetric(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint m = A->m ;
    NLuint i,ij ;
    NLCoeff* c = NULL ;
    for(i=0; i<m; i++) {
        NLRowColumn* Ri = &(A->row[i]) ;
        y[i] = 0 ;
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
            if(i != c->index) {
                y[c->index] += c->value * x[i] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_rows(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    /* 
     * Note: OpenMP does not like unsigned ints
     * (causes some floating point exceptions),
     * therefore I use here signed ints for all
     * indices.
     */
    
    int m = (int)(A->m) ;
    int i,ij ;
    NLCoeff* c = NULL ;
    NLRowColumn* Ri = NULL;
#pragma omp parallel for private(i,ij,c,Ri) 
    for(i=0; i<m; i++) {
        Ri = &(A->row[i]) ;       
        y[i] = 0 ;
        for(ij=0; ij<(int)(Ri->size); ij++) {
            c = &(Ri->coeff[ij]) ;
            y[i] += c->value * x[c->index] ;
        }
    }
}

static void nlSparseMatrix_mult_cols_symmetric(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n ;
    NLuint j,ii ;
    NLCoeff* c = NULL ;
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]) ;       
        y[j] = 0 ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
            if(j != c->index) {
                y[j] += c->value * x[c->index] ;
            }
        }
    }
}

static void nlSparseMatrix_mult_cols(
        NLSparseMatrix* A,
        const NLdouble* x,
        NLdouble* y
) {
    NLuint n = A->n ;
    NLuint j,ii ; 
    NLCoeff* c = NULL ;
    NL_CLEAR_ARRAY(NLdouble, y, A->m) ;
    for(j=0; j<n; j++) {
        NLRowColumn* Cj = &(A->column[j]) ;
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]) ;
            y[c->index] += c->value * x[j] ;
        }
    }
}


/* SparseMatrix x Vector routines, main driver routine */

void nlSparseMatrixMult(NLSparseMatrix* A, const NLdouble* x, NLdouble* y) {
    if(A->storage & NL_MATRIX_STORE_COMPRESSED) {
        nlCRSMatrixMult(A->compressed,x,y);
    } else if(A->storage & NL_MATRIX_STORE_ROWS) {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_rows_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_rows(A, x, y) ;
        }
    } else {
        if(A->storage & NL_MATRIX_STORE_SYMMETRIC) {
            nlSparseMatrix_mult_cols_symmetric(A, x, y) ;
        } else {
            nlSparseMatrix_mult_cols(A, x, y) ;
        }
    }
}


/******* extracted from nl_context.c *******/


NLContextStruct* nlCurrentContext = NULL ;

void nlMatrixVectorProd_default(const NLdouble* x, NLdouble* y) {
    NLSparseMatrix* M = &(nlCurrentContext->M);
    nlSparseMatrixMult(M, x, y) ;
    nlCurrentContext->flops += (NLulong)(nlSparseMatrixNNZ(M)*2);
}

NLContext nlNewContext() {
    NLContextStruct* result     = NL_NEW(NLContextStruct) ;
    result->state               = NL_STATE_INITIAL ;
    result->solver              = NL_SOLVER_DEFAULT ;
    result->max_iterations      = 100 ;
    result->threshold           = 1e-6 ;
    result->omega               = 1.5 ;
    result->row_scaling         = 1.0 ;
    result->right_hand_side     = 0.0 ;
    result->right_hand_side_set = NL_FALSE ;
    result->inner_iterations    = 5 ;
    result->matrix_vector_prod  = nlMatrixVectorProd_default ;
    result->solver_func         = nlDefaultSolver ;
    result->progress_func       = NULL;
    result->verbose             = NL_FALSE;
    nlMakeCurrent(result) ;
    return result ;
}

void nlDeleteContext(NLContext context_in) {
    NLContextStruct* context = (NLContextStruct*)(context_in) ;
    if(nlCurrentContext == context) {
        nlCurrentContext = NULL ;
    }
    if(context->alloc_M) {
        nlSparseMatrixDestroy(&context->M) ;
    }
    if(context->alloc_af) {
        nlRowColumnDestroy(&context->af) ;
    }
    if(context->alloc_al) {
        nlRowColumnDestroy(&context->al) ;
    }
    if(context->alloc_xl) {
        nlRowColumnDestroy(&context->xl) ;
    }
    if(context->alloc_variable) {
        NL_DELETE_ARRAY(context->variable) ;
    }
    if(context->alloc_x) {
        NL_DELETE_ARRAY(context->x) ;
    }
    if(context->alloc_b) {
        NL_DELETE_ARRAY(context->b) ;
    }

#ifdef NL_PARANOID
    NL_CLEAR(NLContextStruct, context) ;
#endif
    NL_DELETE(context) ;
}

void nlMakeCurrent(NLContext context) {
    nlCurrentContext = (NLContextStruct*)(context) ;
}

NLContext nlGetCurrent() {
    return nlCurrentContext ;
}


/* Finite state automaton   */

void nlCheckState(NLenum state) {
    nl_assert(nlCurrentContext->state == state) ;
}

void nlTransition(NLenum from_state, NLenum to_state) {
    nlCheckState(from_state) ;
    nlCurrentContext->state = to_state ;
}


/* Preconditioner setup and default solver */

static void nlSetupPreconditioner() {
    switch(nlCurrentContext->preconditioner) {
    case NL_PRECOND_NONE:
        nlCurrentContext->precond_vector_prod = NULL ;
        break ;
    case NL_PRECOND_JACOBI:
        nlCurrentContext->precond_vector_prod = nlPreconditioner_Jacobi ;
        break ;
    case NL_PRECOND_SSOR:
        nlCurrentContext->precond_vector_prod = nlPreconditioner_SSOR ;
        break ;
    default:
        nl_assert_not_reached ;
        break ;
    }
    /* Check compatibility between solver and preconditioner */
    if(
        nlCurrentContext->solver == NL_BICGSTAB && 
        nlCurrentContext->preconditioner == NL_PRECOND_SSOR
    ) {
        nlWarning(
            "nlSolve", 
            "cannot use SSOR preconditioner with non-symmetric matrix, switching to Jacobi"
        ) ;
        nlCurrentContext->preconditioner = NL_PRECOND_JACOBI ;        
        nlCurrentContext->precond_vector_prod = nlPreconditioner_Jacobi ;
    }
    if(
        nlCurrentContext->solver == NL_GMRES && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for GMRES") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for SUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_PERM_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for PERMSUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
    if(
        nlCurrentContext->solver == NL_SYMMETRIC_SUPERLU_EXT && 
        nlCurrentContext->preconditioner != NL_PRECOND_NONE
    ) {
        nlWarning("nlSolve", "Preconditioner not implemented yet for PERMSUPERLU") ;
        nlCurrentContext->preconditioner = NL_PRECOND_NONE ;        
        nlCurrentContext->precond_vector_prod = NULL ;
    }
}

NLboolean nlDefaultSolver() {
    NLboolean result = NL_TRUE ;
    nlSetupPreconditioner() ;
    switch(nlCurrentContext->solver) {
    case NL_CG: {
        if(nlCurrentContext->preconditioner == NL_PRECOND_NONE) {
            nlCurrentContext->used_iterations = nlSolve_CG() ;
        } else {
            nlCurrentContext->used_iterations = nlSolve_CG_precond() ;
        }
    } break ;
    case NL_BICGSTAB: {
        if(nlCurrentContext->preconditioner == NL_PRECOND_NONE) {
            nlCurrentContext->used_iterations = nlSolve_BICGSTAB() ;
        } else {
            nlCurrentContext->used_iterations = nlSolve_BICGSTAB_precond() ;
        }
    } break ;
    case NL_GMRES: {
        nlCurrentContext->used_iterations = nlSolve_GMRES() ;
    } break ;
    case NL_CNC_FLOAT_CRS_EXT:
    case NL_CNC_DOUBLE_CRS_EXT:
    case NL_CNC_FLOAT_BCRS2_EXT:
    case NL_CNC_DOUBLE_BCRS2_EXT:
    case NL_CNC_FLOAT_ELL_EXT:
    case NL_CNC_DOUBLE_ELL_EXT:
    case NL_CNC_FLOAT_HYB_EXT:
    case NL_CNC_DOUBLE_HYB_EXT: {
        nlCurrentContext->used_iterations = nlSolve_CNC() ;
    } break ;
    case NL_SUPERLU_EXT: 
    case NL_PERM_SUPERLU_EXT: 
    case NL_SYMMETRIC_SUPERLU_EXT: {
        result = nlSolve_SUPERLU() ;
    } break ;
    default:
        nl_assert_not_reached ;
    }
    return result ;
}

/******* extracted from nl_blas.c *******/


#ifdef NL_USE_ATLAS
int NL_FORTRAN_WRAP(xerbla)(char *srname, int *info) {
    printf("** On entry to %6s, parameter number %2d had an illegal value\n",
              srname, *info
    );
    return 0;
} 
#ifndef NL_USE_BLAS
#define NL_USE_BLAS
#endif
#endif

#ifdef NL_USE_SUPERLU
#ifndef NL_USE_BLAS
#define NL_USE_BLAS
/* 
 * The BLAS included in SuperLU does not have DTPSV,
 * we use the DTPSV embedded in OpenNL.
 */
#define NEEDS_DTPSV
#endif
#endif

#ifndef NL_USE_BLAS
#define NEEDS_DTPSV
#endif



/* BLAS routines                                                           */
/* copy-pasted from CBLAS (i.e. generated from f2c) */

/*
 * lsame
 * xerbla
 * daxpy
 * ddot
 * dscal
 * dnrm2
 * dcopy
 * dgemv
 * dtpsv
 */



typedef NLint     integer ;
typedef NLdouble  doublereal ;
typedef NLboolean logical ;
typedef NLint     ftnlen ;


#ifndef max
#define max(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef NL_USE_BLAS

int NL_FORTRAN_WRAP(lsame)(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   

    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
*/  

    /* System generated locals */
    int ret_val;
    
    /* Local variables */
    int inta, intb, zcode;

    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
        return ret_val;
    }

    /* Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

    /* Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {
        /* ASCII is assumed - ZCODE is the ASCII code of either lower or   
          upper case 'Z'. */
        if (inta >= 97 && inta <= 122) inta += -32;
        if (intb >= 97 && intb <= 122) intb += -32;

    } else if (zcode == 233 || zcode == 169) {
        /* EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or   
          upper case 'Z'. */
        if ((inta >= 129 && inta <= 137) || 
            (inta >= 145 && inta <= 153) || 
            (inta >= 162 && inta <= 169)
        )
            inta += 64;
        if (
            (intb >= 129 && intb <= 137) || 
            (intb >= 145 && intb <= 153) || 
            (intb >= 162 && intb <= 169)
        )
            intb += 64;
    } else if (zcode == 218 || zcode == 250) {
        /* ASCII is assumed, on Prime machines - ZCODE is the ASCII code   
          plus 128 of either lower or upper case 'Z'. */
        if (inta >= 225 && inta <= 250) inta += -32;
        if (intb >= 225 && intb <= 250) intb += -32;
    }
    ret_val = inta == intb;
    return ret_val;
    
} /* lsame_ */

/* Subroutine */ int NL_FORTRAN_WRAP(xerbla)(char *srname, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INT   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2d had an illegal value\n",
                srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */


/* Subroutine */ int NL_FORTRAN_WRAP(daxpy)(integer *n, doublereal *da, doublereal *dx, 
        integer *incx, doublereal *dy, integer *incy)
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*da == 0.) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        DY(iy) += *da * DX(ix);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
        DY(i) += *da * DX(i);
        DY(i + 1) += *da * DX(i + 1);
        DY(i + 2) += *da * DX(i + 2);
        DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    nl_arg_used(i__1);
    return 0;
} /* daxpy_ */
#undef DY
#undef DX


doublereal NL_FORTRAN_WRAP(ddot)(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
        integer *incy)
{

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]

    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        dtemp += DX(ix) * DY(iy);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        dtemp += DX(i) * DY(i);
/* L30: */
    }
    if (*n < 5) {
        goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 5) {
        dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) * 
                DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
/* L50: */
    }
L60:
    ret_val = dtemp;
    nl_arg_used(i__1);
    return ret_val;
} /* ddot_ */
#undef DY
#undef DX

/* Subroutine */ int NL_FORTRAN_WRAP(dscal)(integer *n, doublereal *da, doublereal *dx, 
    integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#ifdef DX
#undef DX
#endif
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
        return 0;
    }
    if (*incx == 1) {
        goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
        DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
        DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
        DX(i) = *da * DX(i);
        DX(i + 1) = *da * DX(i + 1);
        DX(i + 2) = *da * DX(i + 2);
        DX(i + 3) = *da * DX(i + 3);
        DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    nl_arg_used(i__1);
    nl_arg_used(i__2);
    return 0;
} /* dscal_ */
#undef DX

doublereal NL_FORTRAN_WRAP(dnrm2)(integer *n, doublereal *x, integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    /* BL: already declared in the included <math.h>, 
       we do not need it here. */
    /*double sqrt(doublereal); */

    /* Local variables */
    static doublereal norm, scale, absxi;
    static integer ix;
    static doublereal ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#ifdef X
#undef X
#endif
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
        norm = 0.;
    } else if (*n == 1) {
        norm = fabs(X(1));
    } else {
        scale = 0.;
        ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
            if (X(ix) != 0.) {
                absxi = (d__1 = X(ix), fabs(d__1));
                if (scale < absxi) {
/* Computing 2nd power */
                    d__1 = scale / absxi;
                    ssq = ssq * (d__1 * d__1) + 1.;
                    scale = absxi;
                } else {
/* Computing 2nd power */
                    d__1 = absxi / scale;
                    ssq += d__1 * d__1;
                }
            }
/* L10: */
        }
        norm = scale * sqrt(ssq);
    }

    ret_val = norm;

    nl_arg_used(i__1);
    nl_arg_used(i__2);

    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */
#undef X

/* Subroutine */ int NL_FORTRAN_WRAP(dcopy)(integer *n, doublereal *dx, integer *incx, 
        doublereal *dy, integer *incy)
{

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
        DY(iy) = DX(ix);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        DY(i) = DX(i);
/* L30: */
    }
    if (*n < 7) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 7) {
        DY(i) = DX(i);
        DY(i + 1) = DX(i + 1);
        DY(i + 2) = DX(i + 2);
        DY(i + 3) = DX(i + 3);
        DY(i + 4) = DX(i + 4);
        DY(i + 5) = DX(i + 5);
        DY(i + 6) = DX(i + 6);
/* L50: */
    }
    nl_arg_used(i__1);
    return 0;
} /* dcopy_ */

#undef DX
#undef DY

/* Subroutine */ int NL_FORTRAN_WRAP(dgemv)(char *trans, integer *m, integer *n, doublereal *
        alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
        doublereal *beta, doublereal *y, integer *incy)
{


    /* System generated locals */
    /* integer a_dim1, a_offset ; */
    integer i__1, i__2; 

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i, j;
/*    extern logical lsame_(char *, char *); */
    static integer ix, iy, jx, jy, kx, ky;
/*    extern int xerbla_(char *, integer *); */


/*  Purpose   
    =======   

    DGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

    Parameters   
    ==========   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! NL_FORTRAN_WRAP(lsame)(trans, "N") && ! NL_FORTRAN_WRAP(lsame)(trans, "T") && ! 
            NL_FORTRAN_WRAP(lsame)(trans, "C")) {
        info = 1;
    } else if (*m < 0) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*lda < max(1,*m)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    } else if (*incy == 0) {
        info = 11;
    }
    if (info != 0) {
        NL_FORTRAN_WRAP(xerbla)("DGEMV ", &info);
        return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
        return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {
        lenx = *n;
        leny = *m;
    } else {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
        if (*incy == 1) {
            if (*beta == 0.) {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(i) = 0.;
/* L10: */
                }
            } else {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(i) = *beta * Y(i);
/* L20: */
                }
            }
        } else {
            iy = ky;
            if (*beta == 0.) {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(iy) = 0.;
                    iy += *incy;
/* L30: */
                }
            } else {
                i__1 = leny;
                for (i = 1; i <= leny; ++i) {
                    Y(iy) = *beta * Y(iy);
                    iy += *incy;
/* L40: */
                }
            }
        }
    }
    if (*alpha == 0.) {
        return 0;
    }
    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (X(jx) != 0.) {
                    temp = *alpha * X(jx);
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        Y(i) += temp * A(i,j);
/* L50: */
                    }
                }
                jx += *incx;
/* L60: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (X(jx) != 0.) {
                    temp = *alpha * X(jx);
                    iy = ky;
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        Y(iy) += temp * A(i,j);
                        iy += *incy;
/* L70: */
                    }
                }
                jx += *incx;
/* L80: */
            }
        }
    } else {

/*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                temp = 0.;
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp += A(i,j) * X(i);
/* L90: */
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
/* L100: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                temp = 0.;
                ix = kx;
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp += A(i,j) * X(ix);
                    ix += *incx;
/* L110: */
                }
                Y(jy) += *alpha * temp;
                jy += *incy;
/* L120: */
            }
        }
    }

    nl_arg_used(i__1);
    nl_arg_used(i__2);
    return 0;

/*     End of DGEMV . */

} /* dgemv_ */

#undef X
#undef Y
#undef A




#else

extern void NL_FORTRAN_WRAP(daxpy)( 
    int *n, double *alpha, double *x,
    int *incx, double *y, int *incy 
) ;

extern double NL_FORTRAN_WRAP(ddot)( 
    int *n, double *x, int *incx, double *y,
    int *incy 
) ;

extern double NL_FORTRAN_WRAP(dnrm2)( int *n, double *x, int *incx ) ;

extern int NL_FORTRAN_WRAP(dcopy)(int* n, double* dx, int* incx, double* dy, int* incy) ;

extern void NL_FORTRAN_WRAP(dscal)(int* n, double* alpha, double *x, int* incx) ;

#ifndef NEEDS_DTPSV
extern void NL_FORTRAN_WRAP(dtpsv)( 
    char *uplo, char *trans, char *diag,
    int *n, double *AP, double *x, int *incx 
) ;
#endif

extern void NL_FORTRAN_WRAP(dgemv)( 
    char *trans, int *m, int *n,
    double *alpha, double *A, int *ldA,
    double *x, int *incx,
    double *beta, double *y, int *incy 
) ;

#endif

#ifdef NEEDS_DTPSV

/* DECK DTPSV */
/* Subroutine */ 
int NL_FORTRAN_WRAP(dtpsv)(
   char* uplo, 
   char* trans, 
   char* diag, 
   integer* n, 
   doublereal* ap, 
   doublereal* x, 
   integer* incx
) {
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
/*    extern logical lsame_(); */
    static integer kk, ix, jx, kx;
/*    extern int xerbla_(); */
    static logical nounit;

/* ***BEGIN PROLOGUE  DTPSV */
/* ***PURPOSE  Solve one of the systems of equations. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B4 */
/* ***TYPE      DOUBLE PRECISION (STPSV-S, DTPSV-D, CTPSV-C) */
/* ***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA */
/* ***AUTHOR  Dongarra, J. J., (ANL) */
/*           Du Croz, J., (NAG) */
/*           Hammarling, S., (NAG) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*  DTPSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix, supplied in packed form. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1))/2). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten */
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/* ***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and */
/*                 Hanson, R. J.  An extended set of Fortran basic linear */
/*                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1, */
/*                 pp. 1-17, March 1988. */
/* ***ROUTINES CALLED  LSAME, XERBLA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   861022  DATE WRITTEN */
/*   910605  Modified to meet SLATEC prologue standards.  Only comment */
/*           lines were modified.  (BKS) */
/* ***END PROLOGUE  DTPSV */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DTPSV */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (!NL_FORTRAN_WRAP(lsame)(uplo, "U") && 
        !NL_FORTRAN_WRAP(lsame)(uplo, "L")
    ) {
        info = 1;
    } else if (
        !NL_FORTRAN_WRAP(lsame)(trans, "N") && 
        !NL_FORTRAN_WRAP(lsame)(trans, "T") && 
        !NL_FORTRAN_WRAP(lsame)(trans, "C")
    ) {
        info = 2;
    } else if (
        !NL_FORTRAN_WRAP(lsame)(diag, "U") && 
        !NL_FORTRAN_WRAP(lsame)(diag, "N")
    ) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*incx == 0) {
        info = 7;
    }
    if (info != 0) {
        NL_FORTRAN_WRAP(xerbla)("DTPSV ", &info);
        return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
        return 0;
    }

    nounit = (logical)(NL_FORTRAN_WRAP(lsame)(diag, "N"));

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (NL_FORTRAN_WRAP(lsame)(trans, "N")) {

/*        Form  x := inv( A )*x. */

        if (NL_FORTRAN_WRAP(lsame)(uplo, "U")) {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.) {
                        if (nounit) {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk - 1;
                        for (i__ = j - 1; i__ >= 1; --i__) {
                            x[i__] -= temp * ap[k];
                            --k;
/* L10: */
                        }
                    }
                    kk -= j;
/* L20: */
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.) {
                        if (nounit) {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__1 = kk - j + 1;
                        for (k = kk - 1; k >= i__1; --k) {
                            ix -= *incx;
                            x[ix] -= temp * ap[k];
/* L30: */
                        }
                    }
                    jx -= *incx;
                    kk -= j;
/* L40: */
                }
            }
        } else {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.) {
                        if (nounit) {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk + 1;
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            x[i__] -= temp * ap[k];
                            ++k;
/* L50: */
                        }
                    }
                    kk += *n - j + 1;
/* L60: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.) {
                        if (nounit) {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__2 = kk + *n - j;
                        for (k = kk + 1; k <= i__2; ++k) {
                            ix += *incx;
                            x[ix] -= temp * ap[k];
/* L70: */
                        }
                    }
                    jx += *incx;
                    kk += *n - j + 1;
/* L80: */
                }
            }
        }
    } else {

/*        Form  x := inv( A' )*x. */

        if (NL_FORTRAN_WRAP(lsame)(uplo, "U")) {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    k = kk;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        temp -= ap[k] * x[i__];
                        ++k;
/* L90: */
                    }
                    if (nounit) {
                        temp /= ap[kk + j - 1];
                    }
                    x[j] = temp;
                    kk += j;
/* L100: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = kx;
                    i__2 = kk + j - 2;
                    for (k = kk; k <= i__2; ++k) {
                        temp -= ap[k] * x[ix];
                        ix += *incx;
/* L110: */
                    }
                    if (nounit) {
                        temp /= ap[kk + j - 1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    kk += j;
/* L120: */
                }
            }
        } else {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    k = kk;
                    i__1 = j + 1;
                    for (i__ = *n; i__ >= i__1; --i__) {
                        temp -= ap[k] * x[i__];
                        --k;
/* L130: */
                    }
                    if (nounit) {
                        temp /= ap[kk - *n + j];
                    }
                    x[j] = temp;
                    kk -= *n - j + 1;
/* L140: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = kx;
                    i__1 = kk - (*n - (j + 1));
                    for (k = kk; k >= i__1; --k) {
                        temp -= ap[k] * x[ix];
                        ix -= *incx;
/* L150: */
                    }
                    if (nounit) {
                        temp /= ap[kk - *n + j];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    kk -= *n - j + 1;
/* L160: */
                }
            }
        }
    }

    return 0;

/*     End of DTPSV . */

} /* dtpsv_ */

#endif 



/* C wrappers for BLAS routines */

/* x <- a*x */
void dscal( int n, double alpha, double *x, int incx ) {
    NL_FORTRAN_WRAP(dscal)(&n,&alpha,x,&incx);
    nlCurrentContext->flops += (NLulong)(n);
}

/* y <- x */
void dcopy( 
    int n, const double *x, int incx, double *y, int incy 
) {
    NL_FORTRAN_WRAP(dcopy)(&n,(double*)x,&incx,y,&incy);
}

/* y <- a*x+y */
void daxpy( 
    int n, double alpha, const double *x, int incx, double *y,
    int incy 
) {
    NL_FORTRAN_WRAP(daxpy)(&n,&alpha,(double*)x,&incx,y,&incy);
    nlCurrentContext->flops += (NLulong)(2*n);    
}

/* returns x^T*y */
double ddot( 
    int n, const double *x, int incx, const double *y, int incy 
) {
    nlCurrentContext->flops += (NLulong)(2*n);        
    return NL_FORTRAN_WRAP(ddot)(&n,(double*)x,&incx,(double*)y,&incy);
}

/* returns |x|_2 */
double dnrm2( int n, const double *x, int incx ) {
    nlCurrentContext->flops += (NLulong)(2*n);        
    return NL_FORTRAN_WRAP(dnrm2)(&n,(double*)x,&incx);
}

/* x <- A^{-1}*x,  x <- A^{-T}*x */
void dtpsv( 
    MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
) {
    static char *UL[2] = { "U", "L" };
    static char *T[3]  = { "N", "T", 0 };
    static char *D[2]  = { "U", "N" };
    /* TODO: update flops */
    NL_FORTRAN_WRAP(dtpsv)(UL[(int)uplo],T[(int)trans],D[(int)diag],&n,(double*)AP,x,&incx); 
}

/* y <- alpha*A*x + beta*y,  y <- alpha*A^T*x + beta*y,   A-(m,n) */
void dgemv( 
    MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
) {
    static char *T[3] = { "N", "T", 0 };
    /* TODO: update flops */    
    NL_FORTRAN_WRAP(dgemv)(T[(int)trans],&m,&n,&alpha,(double*)A,&ldA,(double*)x,&incx,&beta,y,&incy);
}


/* End of BLAS routines */


/******* extracted from nl_iterative_solvers.c *******/



/* Solvers */

/*
 * The implementation of the solvers is inspired by 
 * the lsolver library, by Christian Badura, available from:
 * http://www.mathematik.uni-freiburg.de
 * /IAM/Research/projectskr/lin_solver/
 *
 * About the Conjugate Gradient, details can be found in:
 *  Ashby, Manteuffel, Saylor
 *     A taxononmy for conjugate gradient methods
 *     SIAM J Numer Anal 27, 1542-1568 (1990)
 */

NLuint nlSolve_CG() {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = (NLint)(nlCurrentContext->n) ;

    NLdouble *g = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *r = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *p = NL_NEW_ARRAY(NLdouble, N) ;
    NLuint its=0;
    NLint i;
    NLdouble t, tau, sig, rho, gam;
    NLdouble b_square=ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble accu =0.0;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble curr_err;
    
    nlCurrentContext->matrix_vector_prod(x,g);
    daxpy(N,-1.,b,1,g,1);
    dscal(N,-1.,g,1);
    dcopy(N,g,1,r,1);
    curr_err = ddot(N,g,1,g,1);
    while ( curr_err >err && its < max_iter) {
        if(nlCurrentContext->progress_func != NULL) {
            nlCurrentContext->progress_func(its, max_iter, curr_err, err);
        }
        if(nlCurrentContext->verbose && !(its % 100)) {
            printf ( "%d : %.10e -- %.10e\n", its, curr_err, err ) ;
        }
        nlCurrentContext->matrix_vector_prod(r,p);
        rho=ddot(N,p,1,p,1);
        sig=ddot(N,r,1,p,1);
        tau=ddot(N,g,1,r,1);
        t=tau/sig;
        daxpy(N,t,r,1,x,1);
        daxpy(N,-t,p,1,g,1);
        gam=(t*t*rho-tau)/tau;
        dscal(N,gam,r,1);
        daxpy(N,1.,g,1,r,1);
        ++its;
        curr_err = ddot(N,g,1,g,1);
    }
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i) {
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }
    if(b_square == 0.0) {
        nlCurrentContext->error = sqrt(accu);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
        }
    } else {
        nlCurrentContext->error = sqrt(accu/b_square);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b||/||b|| = %e\n",nlCurrentContext->error);
        }
    }
    
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(g) ;
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(p) ;

    nlCurrentContext->used_iterations = its;
    
    return its;
} 


NLuint nlSolve_CG_precond()  {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = (NLint)(nlCurrentContext->n) ;

    NLdouble* r = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble* d = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble* h = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *Ad = h;
    NLuint its=0;
    NLdouble rh, alpha, beta;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLint i;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;
    NLdouble curr_err;
    
    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    nlCurrentContext->precond_vector_prod(r,d);
    dcopy(N,d,1,h,1);
    rh=ddot(N,r,1,h,1);
    curr_err = ddot(N,r,1,r,1);

    while ( curr_err >err && its < max_iter) {
        if(nlCurrentContext->progress_func != NULL) {
            nlCurrentContext->progress_func(its, max_iter, curr_err, err);
        }
        if( nlCurrentContext->verbose && !(its % 100)) {
           printf ( "%d : %.10e -- %.10e\n", its, curr_err, err ) ;
        }
        nlCurrentContext->matrix_vector_prod(d,Ad);
        alpha=rh/ddot(N,d,1,Ad,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-alpha,Ad,1,r,1);
        nlCurrentContext->precond_vector_prod(r,h);
        beta=1./rh; rh=ddot(N,r,1,h,1); beta*=rh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        ++its;
        curr_err = ddot(N,r,1,r,1);
    }
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i) {
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }

    if(b_square == 0.0) {
        nlCurrentContext->error = sqrt(accu);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
        }
    } else {
        nlCurrentContext->error = sqrt(accu/b_square);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b||/||b|| = %e\n",nlCurrentContext->error);
        }
    }
    
    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(d) ;
    NL_DELETE_ARRAY(h) ;

    nlCurrentContext->used_iterations = its;
    
    return its;
}

NLuint nlSolve_BICGSTAB() {
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = (NLint)(nlCurrentContext->n) ;
    NLint     i ;

    NLdouble *rT  = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *d   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *h   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *u   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *Ad  = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *t   = NL_NEW_ARRAY(NLdouble, N) ; 
    NLdouble *s   = h;
    NLdouble rTh, rTAd, rTr, alpha, beta, omega, st, tt ;
    NLuint its=0;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err=eps*eps*b_square;
    NLdouble *r = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;

    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    dcopy(N,r,1,d,1);
    dcopy(N,d,1,h,1);
    dcopy(N,h,1,rT,1);
    nl_assert( ddot(N,rT,1,rT,1)>1e-40 );
    rTh=ddot(N,rT,1,h,1);
    rTr=ddot(N,r,1,r,1);

    while ( rTr>err && its < max_iter) {
        if(nlCurrentContext->progress_func != NULL) {
            nlCurrentContext->progress_func(its, max_iter, rTr, err);
        }
        if( (nlCurrentContext->verbose) && !(its % 100)) {
            printf ( "%d : %.10e -- %.10e\n", its, rTr, err ) ;
        }
        nlCurrentContext->matrix_vector_prod(d,Ad);
        rTAd=ddot(N,rT,1,Ad,1);
        nl_assert( fabs(rTAd)>1e-40 );
        alpha=rTh/rTAd;
        daxpy(N,-alpha,Ad,1,r,1);
        dcopy(N,h,1,s,1);
        daxpy(N,-alpha,Ad,1,s,1);
        nlCurrentContext->matrix_vector_prod(s,t);
        daxpy(N,1.,t,1,u,1);
        dscal(N,alpha,u,1);
        st=ddot(N,s,1,t,1);
        tt=ddot(N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        daxpy(N,-omega,t,1,r,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-omega,s,1,x,1);
        dcopy(N,s,1,h,1);
        daxpy(N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        daxpy(N,-beta*omega,Ad,1,d,1);
        rTr=ddot(N,r,1,r,1);
        ++its;
    }

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i){
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }

    if(b_square == 0.0) {
        nlCurrentContext->error = sqrt(accu);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
        }
    } else {
        nlCurrentContext->error = sqrt(accu/b_square);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b||/||b|| = %e\n",nlCurrentContext->error);
        }
    }

    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(rT) ;
    NL_DELETE_ARRAY(d) ;
    NL_DELETE_ARRAY(h) ;
    NL_DELETE_ARRAY(u) ;
    NL_DELETE_ARRAY(Ad) ;
    NL_DELETE_ARRAY(t) ;

    nlCurrentContext->used_iterations = its;
    
    return its;
}


NLuint nlSolve_BICGSTAB_precond() {

    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLint     N        = (NLint)(nlCurrentContext->n) ;
    NLint     i;

    NLdouble *rT  = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *d   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *h   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *u   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *Sd  = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *t   = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *aux = NL_NEW_ARRAY(NLdouble, N) ;
    NLdouble *s   = h;
    NLdouble rTh, rTSd, rTr, alpha, beta, omega, st, tt;
    NLuint its=0;
    NLdouble b_square = ddot(N,b,1,b,1);
    NLdouble err  = eps*eps*b_square;
    NLdouble *r   = NL_NEW_ARRAY(NLdouble, N);
    NLdouble * Ax = NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;

    nlCurrentContext->matrix_vector_prod(x,r);
    daxpy(N,-1.,b,1,r,1);
    nlCurrentContext->precond_vector_prod(r,d);
    dcopy(N,d,1,h,1);
    dcopy(N,h,1,rT,1);
    nl_assert( ddot(N,rT,1,rT,1)>1e-40 );
    rTh=ddot(N,rT,1,h,1);
    rTr=ddot(N,r,1,r,1);

    while ( rTr>err && its < max_iter) {
        if(nlCurrentContext->progress_func != NULL) {
            nlCurrentContext->progress_func(its, max_iter, rTr, err);
        }
        if( (nlCurrentContext->verbose) && !(its % 100)) {
            printf ( "%d : %.10e -- %.10e\n", its, rTr, err ) ;
        }
        nlCurrentContext->matrix_vector_prod(d,aux);
        nlCurrentContext->precond_vector_prod(aux,Sd);
        rTSd=ddot(N,rT,1,Sd,1);
        nl_assert( fabs(rTSd)>1e-40 );
        alpha=rTh/rTSd;
        daxpy(N,-alpha,aux,1,r,1);
        dcopy(N,h,1,s,1);
        daxpy(N,-alpha,Sd,1,s,1);
        nlCurrentContext->matrix_vector_prod(s,aux);
        nlCurrentContext->precond_vector_prod(aux,t);
        daxpy(N,1.,t,1,u,1);
        dscal(N,alpha,u,1);
        st=ddot(N,s,1,t,1);
        tt=ddot(N,t,1,t,1);
        if ( fabs(st)<1e-40 || fabs(tt)<1e-40 ) {
            omega = 0.;
        } else {
            omega = st/tt;
        }
        daxpy(N,-omega,aux,1,r,1);
        daxpy(N,-alpha,d,1,x,1);
        daxpy(N,-omega,s,1,x,1);
        dcopy(N,s,1,h,1);
        daxpy(N,-omega,t,1,h,1);
        beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
        dscal(N,beta,d,1);
        daxpy(N,1.,h,1,d,1);
        daxpy(N,-beta*omega,Sd,1,d,1);
        rTr=ddot(N,r,1,r,1);
        ++its;
    }

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < N ; ++i){
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }

    if(b_square == 0.0) {
        nlCurrentContext->error = sqrt(accu);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
        }
    } else {
        nlCurrentContext->error = sqrt(accu/b_square);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b||/||b|| = %e\n",nlCurrentContext->error);
        }
    }

    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(r);
    NL_DELETE_ARRAY(rT);
    NL_DELETE_ARRAY(d);
    NL_DELETE_ARRAY(h);
    NL_DELETE_ARRAY(u);
    NL_DELETE_ARRAY(Sd);
    NL_DELETE_ARRAY(t);
    NL_DELETE_ARRAY(aux);

    nlCurrentContext->used_iterations = its;
        
    return its;
}

NLuint nlSolve_GMRES() {

    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLint    max_iter  = (NLint)(nlCurrentContext->max_iterations) ;
    NLint    n         = (NLint)(nlCurrentContext->n) ;
    NLint    m         = (NLint)(nlCurrentContext->inner_iterations) ;

    typedef NLdouble *NLdoubleP;
    NLdouble *V   = NL_NEW_ARRAY(NLdouble, n*(m+1)   ) ;
    NLdouble *U   = NL_NEW_ARRAY(NLdouble, m*(m+1)/2 ) ;
    NLdouble *r   = NL_NEW_ARRAY(NLdouble, n         ) ;
    NLdouble *y   = NL_NEW_ARRAY(NLdouble, m+1       ) ;
    NLdouble *c   = NL_NEW_ARRAY(NLdouble, m         ) ;
    NLdouble *s   = NL_NEW_ARRAY(NLdouble, m         ) ;
    NLdouble **v  = NL_NEW_ARRAY(NLdoubleP, m+1      ) ;
    NLdouble * Ax = NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu =0.0;
    NLint i, j, io, uij, u0j ; 
    NLint its = -1 ;
    NLdouble beta, h, rd, dd, nrm2b ;

    for ( i=0; i<=m; ++i ){
        v[i]=V+i*n ;
    }
    
    nrm2b=dnrm2(n,b,1);
    io=0;

    do  { /* outer loop */
        ++io;
        nlCurrentContext->matrix_vector_prod(x,r);
        daxpy(n,-1.,b,1,r,1);
        beta=dnrm2(n,r,1);
        dcopy(n,r,1,v[0],1);
        dscal(n,1./beta,v[0],1);

        y[0]=beta;
        j=0;
        uij=0;
        do { /* inner loop: j=0,...,m-1 */
            u0j=uij;
            nlCurrentContext->matrix_vector_prod(v[j],v[j+1]);
            dgemv(
                Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1
            );
            dgemv(
                NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1
            );
            h=dnrm2(n,v[j+1],1);
            dscal(n,1./h,v[j+1],1);
            for (i=0; i<j; ++i ) { /* rotiere neue Spalte */
                double tmp = c[i]*U[uij]-s[i]*U[uij+1];
                U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
                U[uij]     = tmp;
                ++uij;
            }
            { /* berechne neue Rotation */
                rd     = U[uij];
                dd     = sqrt(rd*rd+h*h);
                c[j]   = rd/dd;
                s[j]   = -h/dd;
                U[uij] = dd;
                ++uij;
            }
            { /* rotiere rechte Seite y (vorher: y[j+1]=0) */
                y[j+1] = s[j]*y[j];
                y[j]   = c[j]*y[j];
            }
            ++j;
        } while ( 
            j<m && fabs(y[j])>=eps*nrm2b 
        ) ;
        { /* minimiere bzgl Y */
            dtpsv(
                UpperTriangle,
                NoTranspose,
                NotUnitTriangular,
                j,U,y,1
            );
            /* correct X */
            dgemv(NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
        }
    } while ( fabs(y[j])>=eps*nrm2b && (m*(io-1)+j) < max_iter);
    
    /* Count the inner iterations */
    its = m*(io-1)+j;

    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < n ; ++i) {
        accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
    }

    if(nrm2b == 0.0) {
        nlCurrentContext->error = sqrt(accu);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b|| = %e\n",nlCurrentContext->error);
        }
    } else {
        nlCurrentContext->error = sqrt(accu/nrm2b);
        if(nlCurrentContext->verbose) {
            printf("in OpenNL : ||Ax-b||/||b|| = %e\n",nlCurrentContext->error);
        }
    }

    NL_DELETE_ARRAY(Ax);
    NL_DELETE_ARRAY(V) ;
    NL_DELETE_ARRAY(U) ;
    NL_DELETE_ARRAY(r) ;
    NL_DELETE_ARRAY(y) ;
    NL_DELETE_ARRAY(c) ;
    NL_DELETE_ARRAY(s) ;
    NL_DELETE_ARRAY(v) ;

    nlCurrentContext->used_iterations = (NLuint)its;
    
    return (NLuint)its;
}






/******* extracted from nl_preconditioners.c *******/



/* preconditioners */

/* Utilities for preconditioners */

void nlMultDiagonal(NLdouble* xy, NLdouble omega) {
    NLuint N = nlCurrentContext->n ;
    NLuint i ;
    NLdouble* diag = nlCurrentContext->M.diag ;
    for(i=0; i<N; i++) {
        xy[i] *= (diag[i] / omega) ;
    }
    nlCurrentContext->flops += (NLulong)(N);
}

void nlMultDiagonalInverse(NLdouble* xy, NLdouble omega) {
    NLuint N = nlCurrentContext->n ;
    NLuint i ;
    NLdouble* diag = nlCurrentContext->M.diag ;
    for(i=0; i<N; i++) {
        xy[i] *= ((diag[i] != 0) ? (omega / diag[i]) : omega) ;
    }
    nlCurrentContext->flops += (NLulong)(N);    
}

void nlMultLowerInverse(const NLdouble* x, NLdouble* y, double omega) {
    NLSparseMatrix* A = &(nlCurrentContext->M) ;
    NLuint n       = A->n ;
    NLdouble* diag = A->diag ;
    NLuint i ;
    NLuint ij ;
    NLCoeff* c = NULL ;
    NLdouble S ;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC) ;
    nl_assert(A->storage & NL_MATRIX_STORE_ROWS) ;

    for(i=0; i<n; i++) {
        NLRowColumn*  Ri = &(A->row[i]) ;       
        S = 0 ;
        for(ij=0; ij < Ri->size; ij++) {
            c = &(Ri->coeff[ij]) ;
            nl_parano_assert(c->index <= i) ; 
            if(c->index != i) {
                S += c->value * y[c->index] ; 
            }
        }
        nlCurrentContext->flops += (NLulong)(2*Ri->size);                    
        y[i] = (x[i] - S) * omega / diag[i] ;
    }
    nlCurrentContext->flops += (NLulong)(n*3);                
}

void nlMultUpperInverse(const NLdouble* x, NLdouble* y, NLdouble omega) {
    NLSparseMatrix* A = &(nlCurrentContext->M) ;
    NLuint n       = A->n ;
    NLdouble* diag = A->diag ;
    NLint i ;
    NLuint ij ;
    NLCoeff* c = NULL ;
    NLdouble S ;

    nl_assert(A->storage & NL_MATRIX_STORE_SYMMETRIC) ;
    nl_assert(A->storage & NL_MATRIX_STORE_COLUMNS) ;

    for(i=(NLint)(n-1); i>=0; i--) {
        NLRowColumn*  Ci = &(A->column[i]) ;       
        S = 0 ;
        for(ij=0; ij < Ci->size; ij++) {
            c = &(Ci->coeff[ij]) ;
            nl_parano_assert(c->index >= i) ; 
            if((NLint)(c->index) != i) {
                S += c->value * y[c->index] ; 
            }
        }
        nlCurrentContext->flops += (NLulong)(2*Ci->size);                    
        y[i] = (x[i] - S) * omega / diag[i] ;
    }
    nlCurrentContext->flops += (NLulong)(n*3);                
}


void nlPreconditioner_Jacobi(const NLdouble* x, NLdouble* y) {
    if(nlCurrentContext->M.storage & NL_MATRIX_STORE_DIAG_INV) {
        NLuint i;
        for(i=0; i<nlCurrentContext->M.diag_size; ++i) {
            y[i] = x[i]*nlCurrentContext->M.diag_inv[i];
        }
        nlCurrentContext->flops += (NLulong)(nlCurrentContext->M.diag_size*3);                        
    } else {
        NLint N = (NLint)(nlCurrentContext->n) ;
        dcopy(N, x, 1, y, 1) ;
        nlMultDiagonalInverse(y, 1.0) ;
    }
}


static double* nlPreconditioner_SSOR_work = NULL;
static NLuint nlPreconditioner_SSOR_work_size = 0;

void nlPreconditioner_SSOR_terminate(void) {
    NL_DELETE_ARRAY(nlPreconditioner_SSOR_work);
}

void nlPreconditioner_SSOR(const NLdouble* x, NLdouble* y) {
    NLdouble omega = nlCurrentContext->omega ;
    NLuint n = nlCurrentContext->n ;
    static NLboolean init = NL_FALSE;
    if(!init) {
        atexit(nlPreconditioner_SSOR_terminate);
        init = NL_TRUE;
    }
    if(n != nlPreconditioner_SSOR_work_size) {
        nlPreconditioner_SSOR_work = NL_RENEW_ARRAY(
            NLdouble, nlPreconditioner_SSOR_work, n
        ) ;
        nlPreconditioner_SSOR_work_size = n ;
    }
    
    nlMultLowerInverse(
        x, nlPreconditioner_SSOR_work, omega
    );
    nlMultDiagonal(
        nlPreconditioner_SSOR_work, omega
    );
    nlMultUpperInverse(
        nlPreconditioner_SSOR_work, y, omega
    );

    dscal((NLint)n, 2.0 - omega, y, 1) ;

    nlCurrentContext->flops += (NLulong)(n);    
}


/******* extracted from nl_superlu.c *******/


/*
 * \file Weak-coupling adapter to call SuperLU from OpenNL, 
 *  works with both SuperLU 3.x and SuperLU 4.x.
 */

#ifdef GEO_DYNAMIC_LIBS
#  ifdef unix
#include <dlfcn.h>

typedef void (*FUNPTR)();

FUNPTR fun_dlsym(void* handle, const char* name) {
    union {
        void* ptr;
        FUNPTR fptr;
    } u;
    u.ptr = dlsym(handle, name);
    return u.fptr;
}
#  endif
#endif








/*
 * Important note: 
 * the order of some constants and the size
 * of some structures have changed between
 * SuperLU 3.x and SuperLU 4.x.
 * See documentation of the SuperLU_version() 
 * function in nl_superlu.c for more details.
 */

typedef enum {
    /* SuperLU 3.x */
    SLU3_NC  =0,   /* column-wise, no supernode */
    SLU3_NR  =1,   /* row-wize, no supernode */
    SLU3_SC  =2,   /* column-wise, supernode */
    SLU3_SR  =3,   /* row-wise, supernode */
    SLU3_NCP =4,   /* column-wise, column-permuted, no supernode*/ 
    SLU3_DN  =5,   /* Fortran style column-wise storage for dense matrix */

    /* SuperLU 4.x */    
    SLU4_NC     =0, /* column-wise, no supernode */
    SLU4_NCP    =1, /* column-wise, column-permuted, no supernode */ 
    SLU4_NR     =2, /* row-wize, no supernode */
    SLU4_SC     =3, /* column-wise, supernode */
    SLU4_SCP    =4, /* supernode, column-wise, permuted */
    SLU4_SR     =5, /* row-wise, supernode */
    SLU4_DN     =6, /* Fortran style column-wise storage for dense matrix */
    SLU4_NR_loc =7  /* distributed compressed row format  */
} Stype_t;

typedef enum {
    SLU_S,     /* single */
    SLU_D,     /* double */
    SLU_C,     /* single complex */
    SLU_Z      /* double complex */
} Dtype_t;

typedef enum {
    SLU_GE,    /* general */
    SLU_TRLU,  /* lower triangular, unit diagonal */
    SLU_TRUU,  /* upper triangular, unit diagonal */
    SLU_TRL,   /* lower triangular */
    SLU_TRU,   /* upper triangular */
    SLU_SYL,   /* symmetric, store lower half */
    SLU_SYU,   /* symmetric, store upper half */
    SLU_HEL,   /* Hermitian, store lower half */
    SLU_HEU    /* Hermitian, store upper half */
} Mtype_t;

typedef int int_t;

typedef struct {
        Stype_t Stype; /* Storage type: interprets the storage structure 
                          pointed to by *Store. */
        Dtype_t Dtype; /* Data type. */
        Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
                          the matrix. */
        int_t  nrow;   /* number of rows */
        int_t  ncol;   /* number of columns */
        void *Store;   /* pointer to the actual storage of the matrix */
} SuperMatrix;

typedef struct {
    int_t lda;    /* leading dimension */
    void *nzval;  /* array of size lda*ncol to represent a dense matrix */
} DNformat;

typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;

typedef enum {
    /* SuperLU 3.x */
    SLU3_NATURAL       = 0,
    SLU3_MMD_ATA       = 1,
    SLU3_MMD_AT_PLUS_A = 2,
    SLU3_COLAMD        = 3,
    SLU3_MY_PERMC      = 4,

    /* SuperLU 4.x */
    SLU4_NATURAL         = 0,
    SLU4_MMD_ATA         = 1,
    SLU4_MMD_AT_PLUS_A   = 2,
    SLU4_COLAMD          = 3,
    SLU4_METIS_AT_PLUS_A = 4,
    SLU4_PARMETIS        = 5,
    SLU4_ZOLTAN          = 6,
    SLU4_MY_PERMC        = 7
} colperm_t;


typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SLU_SINGLE=1, SLU_DOUBLE, SLU_EXTRA}    IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB, SLU4_LLVL, SLU4_ULVL}    MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;
typedef enum {SLU4_ONE_NORM, SLU4_TWO_NORM, SLU4_INF_NORM}      norm_t;
typedef enum {
    SLU4_SILU, SLU4_SMILU_1, SLU4_SMILU_2, SLU4_SMILU_3
} milu_t;

typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    yes_no_t      PrintStat;
    yes_no_t      SymmetricMode;
    double        DiagPivotThresh;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    yes_no_t      ReplaceTinyPivot;
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
} superlu3_options_t;

typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int           ILU_DropRule;
    double        ILU_DropTol;    /* threshold for dropping */
    double        ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t        ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double        ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t        ILU_MILU;
    double        ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    int           nnzL, nnzU;      /* used to store nnzs for now       */
    int           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
                                      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu4_options_t;

typedef void* superlu_options_ptr;

typedef float    flops_t;
typedef unsigned char Logical;

typedef struct {
    int     *panel_histo;    /* histogram of panel size distribution */
    double  *utime;          /* running time at various phases */
    flops_t *ops;            /* operation count at various phases */
    int     TinyPivots;      /* number of tiny pivots */
    int     RefineSteps;     /* number of iterative refinement steps */
    int     slu4_expansions; /* number of memory expansions (SuperLU4) */
} SuperLUStat_t;












typedef void (*FUNPTR_set_default_options)(superlu_options_ptr options);

typedef void (*FUNPTR_StatInit)(SuperLUStat_t *);
typedef void (*FUNPTR_StatFree)(SuperLUStat_t *);

typedef void (*FUNPTR_dCreate_CompCol_Matrix)(
    SuperMatrix *, int, int, int, double *,
    int *, int *, Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_dCreate_Dense_Matrix)(
    SuperMatrix *, int, int, double *, int,
    Stype_t, Dtype_t, Mtype_t);

typedef void (*FUNPTR_Destroy_SuperNode_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_CompCol_Matrix)(SuperMatrix *);
typedef void (*FUNPTR_Destroy_SuperMatrix_Store)(SuperMatrix *);

typedef void (*FUNPTR_dgssv)(
    superlu_options_ptr, SuperMatrix *, int *, int *, SuperMatrix *,
    SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *
);

typedef struct {
    FUNPTR_set_default_options set_default_options;
    FUNPTR_StatInit StatInit;
    FUNPTR_StatFree StatFree;
    FUNPTR_dCreate_CompCol_Matrix dCreate_CompCol_Matrix;
    FUNPTR_dCreate_Dense_Matrix dCreate_Dense_Matrix;
    FUNPTR_Destroy_SuperNode_Matrix Destroy_SuperNode_Matrix;
    FUNPTR_Destroy_CompCol_Matrix Destroy_CompCol_Matrix;
    FUNPTR_Destroy_SuperMatrix_Store Destroy_SuperMatrix_Store;
    FUNPTR_dgssv dgssv;

    void* DLL_handle;

    double version;
} SuperLUContext;

SuperLUContext* SuperLU() {
    static SuperLUContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

/*
 * \brief Gets the version of SuperLU that was dynamically loaded.
 * \details It is important to know the version of SuperLU because
 *  the order of some constants and the size of some structures 
 *  have changed between SuperLU 3.x and SuperLU 4.x.
 *  When there is a mismatch between both versions, constants are 
 *  prefixed by SLU3_ or SLU4_ according to the version.
 * This concerns:
 * - enum constants in Stype_t
 * - enum constants in colperm_t
 * - enum constants in norm_t and milu_t (that only exist in version 4.x)
 * - struct superlu_options_t (use instead superlu3_options_t 
 *  or superlu4_options_t according to version)
 */
double SuperLU_version() {
    return SuperLU()->version;
}

NLboolean SuperLU_is_initialized() {
    return
        SuperLU()->DLL_handle != NULL &&
        SuperLU()->set_default_options != NULL &&
        SuperLU()->StatInit != NULL &&
        SuperLU()->StatFree != NULL &&
        SuperLU()->dCreate_CompCol_Matrix != NULL &&
        SuperLU()->dCreate_Dense_Matrix != NULL &&
        SuperLU()->Destroy_SuperNode_Matrix != NULL &&
        SuperLU()->Destroy_CompCol_Matrix != NULL &&
        SuperLU()->Destroy_SuperMatrix_Store != NULL &&
        SuperLU()->dgssv != NULL ;
}

#define find_superlu_func(name)                                   \
    if(                                                           \
        (                                                         \
            SuperLU()->name =                                     \
            (FUNPTR_##name)fun_dlsym(SuperLU()->DLL_handle,#name) \
        ) == NULL                                                 \
    ) {                                                           \
        nlError("nlInitExtension_SUPERLU/dlsym",dlerror());   \
        return NL_FALSE;                                          \
    }

#define SLU_NR ((SuperLU()->version >= 4.0) ? SLU4_NR : SLU3_NR)

#define SLU_DN ((SuperLU()->version >= 4.0) ? SLU4_DN : SLU3_DN)



NLboolean nlSolve_SUPERLU() {

    /* OpenNL Context */
    NLSparseMatrix* M  = &(nlCurrentContext->M) ;
    NLdouble* b          = nlCurrentContext->b ;
    NLdouble* x          = nlCurrentContext->x ;

    /* Compressed Row Storage matrix representation */
    NLuint    n      = nlCurrentContext->n ;
    NLuint    nnz    = nlSparseMatrixNNZ(M) ; /* Number of Non-Zero coeffs */
    NLint*    xa     = NL_NEW_ARRAY(NLint, n+1) ;
    NLdouble* rhs    = NL_NEW_ARRAY(NLdouble, n) ;
    NLdouble* a      = NL_NEW_ARRAY(NLdouble, nnz) ;
    NLint*    asub   = NL_NEW_ARRAY(NLint, nnz) ;

    /* Permutation vector */
    NLint*    perm_r  = NL_NEW_ARRAY(NLint, n) ;
    NLint*    perm    = NL_NEW_ARRAY(NLint, n) ;

    /* SuperLU variables */
    SuperMatrix A, B ; /* System       */
    SuperMatrix L, U ; /* Factorization of A */
    NLint info ;       /* status code  */
    DNformat *vals = NULL ; /* access to result */
    double *rvals  = NULL ; /* access to result */

    /* SuperLU options and stats */
    superlu3_options_t options3 ;
    superlu4_options_t options4 ;
    SuperLUStat_t     stat ;

    /* Temporary variables */
    NLRowColumn* Ri = NULL ;
    NLuint         i,jj,count ;
    
    /* Sanity checks */
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC)) ;
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS) ;
    nl_assert(M->m == M->n) ;

    if(!SuperLU_is_initialized()) {
        nlError(
            "nlSolve_SUPERLU",
            "SuperLU extension not initialized (nlInitExtension(\"SUPERLU\") missing or failed)"
        );
        return NL_FALSE;
    }
    
    /*
     * Step 1: convert matrix M into SuperLU compressed column 
     *   representation.
     * -------------------------------------------------------
     */

    count = 0 ;
    for(i=0; i<n; i++) {
        Ri = &(M->row[i]) ;
        xa[i] = (NLint)(count);
        for(jj=0; jj<Ri->size; jj++) {
            a[count]    = Ri->coeff[jj].value ;
            asub[count] = (NLint)(Ri->coeff[jj].index) ;
            count++ ;
        }
    }
    xa[n] = (NLint)(nnz) ;

    /* Save memory for SuperLU */
    nlSparseMatrixClear(M) ;


    /*
     * Rem: SuperLU does not support symmetric storage.
     * In fact, for symmetric matrix, what we need 
     * is a SuperLLt algorithm (SuperNodal sparse Cholesky),
     * but it does not exist, anybody wants to implement it ?
     * However, this is not a big problem (SuperLU is just
     * a superset of what we really need.
     */
    SuperLU()->dCreate_CompCol_Matrix(
        &A, (int)n, (int)n, (int)nnz,
        a, asub, xa, 
        SLU_NR,              /* Row_wise, no supernode */
        SLU_D,               /* doubles                */ 
        SLU_GE               /* general storage        */
    );

    /* Step 2: create vector */
    SuperLU()->dCreate_Dense_Matrix(
        &B, (int)n, 1, b, (int)n, 
        SLU_DN, /* Fortran-type column-wise storage */
        SLU_D,  /* doubles                          */
        SLU_GE  /* general                          */
    );
            

    /* Step 3: set SuperLU options 
     * ------------------------------
     */

    if(SuperLU_version() >= 4.0) {
        SuperLU()->set_default_options(&options4) ;
        switch(nlCurrentContext->solver) {
        case NL_SUPERLU_EXT: {
            options4.ColPerm = SLU4_NATURAL ;
        } break ;
        case NL_PERM_SUPERLU_EXT: {
            options4.ColPerm = SLU4_COLAMD ;
        } break ;
        case NL_SYMMETRIC_SUPERLU_EXT: {
            options4.ColPerm = SLU4_MMD_AT_PLUS_A ;
            options4.SymmetricMode = YES ;
        } break ;
        default: {
            nl_assert_not_reached ;
        } break ;
        }
    } else {
        SuperLU()->set_default_options(&options3) ;
        switch(nlCurrentContext->solver) {
        case NL_SUPERLU_EXT: {
            options3.ColPerm = SLU3_NATURAL ;
        } break ;
        case NL_PERM_SUPERLU_EXT: {
            options3.ColPerm = SLU3_COLAMD ;
        } break ;
        case NL_SYMMETRIC_SUPERLU_EXT: {
            options3.ColPerm = SLU3_MMD_AT_PLUS_A ;
            options3.SymmetricMode = YES ;
        } break ;
        default: {
            nl_assert_not_reached ;
        } break ;
        }
    }
    
    SuperLU()->StatInit(&stat) ;

    /* Step 4: call SuperLU main routine
     * ---------------------------------
     */

    if(SuperLU_version() >= 4.0) {
        SuperLU()->dgssv(
            &options4, &A, perm, perm_r, &L, &U, &B, &stat, &info
        );
    } else {
        SuperLU()->dgssv(
            &options3, &A, perm, perm_r, &L, &U, &B, &stat, &info
        );
    }

    /* Step 5: get the solution
     * ------------------------
     * Fortran-type column-wise storage
     */
    vals = (DNformat*)B.Store;
    rvals = (double*)(vals->nzval);
    if(info == 0) {
        for(i = 0; i <  n; i++){
            x[i] = rvals[i];
        }
    } else {
        nlError("nlSolve", "SuperLU failed") ;
    }

    /* Step 6: cleanup
     * ---------------
     */

    /*
     *  For these two ones, only the "store" structure
     * needs to be deallocated (the arrays have been allocated
     * by us).
     */
    SuperLU()->Destroy_SuperMatrix_Store(&A) ;
    SuperLU()->Destroy_SuperMatrix_Store(&B) ;

    /*
     *   These ones need to be fully deallocated (they have been
     * allocated by SuperLU).
     */
    SuperLU()->Destroy_SuperNode_Matrix(&L);
    SuperLU()->Destroy_CompCol_Matrix(&U);

    /* There are some dynamically allocated vectors in the stats */
    SuperLU()->StatFree(&stat) ;

    NL_DELETE_ARRAY(xa) ;
    NL_DELETE_ARRAY(rhs) ;
    NL_DELETE_ARRAY(a) ;
    NL_DELETE_ARRAY(asub) ;
    NL_DELETE_ARRAY(perm_r) ;
    NL_DELETE_ARRAY(perm) ;

    return (info == 0) ;
}


#if defined(GEO_DYNAMIC_LIBS) && defined(unix)

void nlTerminateExtension_SUPERLU() {
    if(SuperLU()->DLL_handle != NULL) {
        dlclose(SuperLU()->DLL_handle);
        SuperLU()->DLL_handle = NULL;
    }
}

NLboolean nlInitExtension_SUPERLU() {
    
    if(SuperLU()->DLL_handle != NULL) {
        return SuperLU_is_initialized();
    }

    SuperLU()->DLL_handle = dlopen("libsuperlu.so", RTLD_NOW);
    if(SuperLU()->DLL_handle == NULL) {
        nlError("nlInitExtension_SUPERLU/dlopen",dlerror());
        return NL_FALSE;
    }

    /* 
     * Check for SuperLU version:
     * Since ILU (incomplete Cholesky) is only available in 4.x, if
     * we find one of the ILU-related symbols in there, then we got
     * a 4.x.
     * TODO: there may be a finer way to detect version.
     */
    if(dlsym(SuperLU()->DLL_handle,"ilu_set_default_options") != NULL) {
        SuperLU()->version = 4.0;
    } else {
        SuperLU()->version = 3.0;
    }
    
    find_superlu_func(set_default_options);
    find_superlu_func(StatInit);
    find_superlu_func(StatFree);
    find_superlu_func(dCreate_CompCol_Matrix);
    find_superlu_func(dCreate_Dense_Matrix);
    find_superlu_func(Destroy_SuperNode_Matrix);
    find_superlu_func(Destroy_CompCol_Matrix);    
    find_superlu_func(Destroy_SuperMatrix_Store);
    find_superlu_func(dgssv);

    atexit(nlTerminateExtension_SUPERLU);
    return NL_TRUE;
}

#else

NLboolean nlInitExtension_SUPERLU() {
    nlError(
        "nlInitExtension_SUPERLU",
        "Only supported in Linux / dynamic lib. mode"
    );
    return NL_FALSE;
}

#endif


/******* extracted from nl_cnc_gpu_cuda.c *******/


NLboolean nlSolverIsCNC(NLint solver){
    return solver == NL_CNC_FLOAT_CRS_EXT 
        || solver == NL_CNC_DOUBLE_CRS_EXT 
        || solver == NL_CNC_FLOAT_BCRS2_EXT 
        || solver == NL_CNC_DOUBLE_BCRS2_EXT         
        || solver == NL_CNC_FLOAT_ELL_EXT         
        || solver == NL_CNC_DOUBLE_ELL_EXT         
        || solver == NL_CNC_FLOAT_HYB_EXT         
        || solver == NL_CNC_DOUBLE_HYB_EXT ;        
}



/* CNC wrapper */

#ifdef NL_USE_CNC

NLuint nlSolve_CNC() {
    unsigned int i;
    NLdouble* b        = nlCurrentContext->b ;
    NLdouble* x        = nlCurrentContext->x ;
    NLdouble  eps      = nlCurrentContext->threshold ;
    NLuint    max_iter = nlCurrentContext->max_iterations ;
    NLSparseMatrix *M  = &(nlCurrentContext->M);
    
    /* local variables for the final error computation */
    NLuint val_ret;
    NLdouble * Ax=NL_NEW_ARRAY(NLdouble,nlCurrentContext->n);
    NLdouble accu     = 0.0;
    NLdouble b_square = 0.0;
    
    
    val_ret=cnc_solve_cg(M, b, x, max_iter, eps, nlCurrentContext->solver);
    
    /* compute the final error */
    nlCurrentContext->matrix_vector_prod(x,Ax);
    for(i = 0 ; i < M->n ; ++i) { 
        accu     +=(Ax[i]-b[i])*(Ax[i]-b[i]);
        b_square += b[i]*b[i]; 
    } 
    printf("in OpenNL : ||Ax-b||/||b|| = %e\n",sqrt(accu)/sqrt(b_square));
    /* cleaning */
    NL_DELETE_ARRAY(Ax);
    return val_ret;
}

#else

NLuint nlSolve_CNC() {
    nl_assert_not_reached ;
    return 0;
}

#endif

/******* extracted from nl_api.c *******/


NLboolean nlInitExtension(const char* extension) {

    nl_arg_used(extension);

    if(!strcmp(extension, "SUPERLU")) {
        return nlInitExtension_SUPERLU() ;
    }

#ifdef NL_USE_CNC
    if(!strcmp(extension, "CNC")) {
        return NL_TRUE ;
    }
#endif
    return NL_FALSE ;
}




/* Get/Set parameters */

void nlSolverParameterd(NLenum pname, NLdouble param) {
    nlCheckState(NL_STATE_INITIAL) ;
    switch(pname) {
    case NL_THRESHOLD: {
        nl_assert(param >= 0) ;
        nlCurrentContext->threshold = (NLdouble)param ;
        nlCurrentContext->threshold_defined = NL_TRUE ;
    } break ;
    case NL_OMEGA: {
        nl_range_assert(param,1.0,2.0) ;
        nlCurrentContext->omega = (NLdouble)param ;
    } break ;
    default: {
        nlError("nlSolverParameterd","Invalid parameter");
        nl_assert_not_reached ;
    } break ;
    }
}

void nlSolverParameteri(NLenum pname, NLint param) {
    nlCheckState(NL_STATE_INITIAL) ;
    switch(pname) {
    case NL_SOLVER: {
        nlCurrentContext->solver = (NLenum)param ;
    } break ;
    case NL_NB_VARIABLES: {
        nl_assert(param > 0) ;
        nlCurrentContext->nb_variables = (NLuint)param ;
    } break ;
    case NL_LEAST_SQUARES: {
        nlCurrentContext->least_squares = (NLboolean)param ;
    } break ;
    case NL_MAX_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->max_iterations = (NLuint)param ;
        nlCurrentContext->max_iterations_defined = NL_TRUE;
    } break ;
    case NL_SYMMETRIC: {
        nlCurrentContext->symmetric = (NLboolean)param ;        
    } break ;
    case NL_INNER_ITERATIONS: {
        nl_assert(param > 0) ;
        nlCurrentContext->inner_iterations = (NLuint)param ;
    } break ;
    case NL_PRECONDITIONER: {
        nlCurrentContext->preconditioner = (NLuint)param ;
        nlCurrentContext->preconditioner_defined = NL_TRUE;
    } break ;
    default: {
        nlError("nlSolverParameteri","Invalid parameter");
        nl_assert_not_reached ;
    } break ;
    }
}


static void nlRowParameterd_obsolete() {
    static NLboolean first_time = NL_TRUE;
    if(first_time) {
        fprintf(stderr,"==== OpenNL warning: called an obsolete function nlRowParameterd() ===========\n");
        fprintf(stderr,"Now use nlRightHandSide() / nlRowScaling() instead\n");
        fprintf(stderr,"   PLEASE TAKE NOTE THAT THE SIGN OF THE ARGUMENT CHANGED IN nlRightHandSide()\n");
        fprintf(stderr,"   AS COMPARED TO nlRowParameterd(NL_RIGHT_HAND_SIDE,.)\n");
        fprintf(stderr,"   (please refer to the documentation and examples in nl.h\n");
        fprintf(stderr,"==============================================================================\n");
        first_time = NL_FALSE;
    }
}

void nlRowParameterd(NLenum pname, NLdouble param) {
    nlRowParameterd_obsolete();    
    nlCheckState(NL_STATE_MATRIX) ;
    switch(pname) {
    case NL_RIGHT_HAND_SIDE: {
        /*
         * Argh: -param, because the old API was inversed,
         * right hand side was b in Ax + b = 0 instead of Ax = b,
         * it is too stupid.
         */
        if(nlCurrentContext->least_squares) {
            nlCurrentContext->right_hand_side =  param ;            
        } else {
            nlCurrentContext->right_hand_side = -param ;
        }
    } break ;
    case NL_ROW_SCALING: {
        nlCurrentContext->row_scaling = param ;
    } break ;
    default: {
        nlError("nlSolverParameterd","Invalid parameter");
        nl_assert_not_reached ;
    }
    }
}

void nlGetBooleanv(NLenum pname, NLboolean* params) {
    switch(pname) {
    case NL_LEAST_SQUARES: {
        *params = nlCurrentContext->least_squares ;
    } break ;
    case NL_SYMMETRIC: {
        *params = nlCurrentContext->symmetric ;
    } break ;
    default: {
        nlError("nlGetBooleanv","Invalid parameter");
        nl_assert_not_reached ;
    } break ;
    }
}

void nlGetDoublev(NLenum pname, NLdouble* params) {
    switch(pname) {
    case NL_THRESHOLD: {
        *params = nlCurrentContext->threshold ;
    } break ;
    case NL_OMEGA: {
        *params = nlCurrentContext->omega ;
    } break ;
    case NL_ERROR: {
        *params = nlCurrentContext->error ;
    } break ;
    case NL_ELAPSED_TIME: {
        *params = nlCurrentContext->elapsed_time ;        
    } break ;
    case NL_GFLOPS: {
        if(nlCurrentContext->elapsed_time == 0) {
            *params = 0.0;
        } else {
            *params = (NLdouble)(nlCurrentContext->flops) /
                (nlCurrentContext->elapsed_time * 1e9) ;
        }
    } break;
    default: {
        nlError("nlGetDoublev","Invalid parameter");
        nl_assert_not_reached ;
    } break ;
    }
}

void nlGetIntegerv(NLenum pname, NLint* params) {
    switch(pname) {
    case NL_SOLVER: {
        *params = (NLint)(nlCurrentContext->solver) ;
    } break ;
    case NL_NB_VARIABLES: {
        *params = (NLint)(nlCurrentContext->nb_variables) ;
    } break ;
    case NL_LEAST_SQUARES: {
        *params = (NLint)(nlCurrentContext->least_squares) ;
    } break ;
    case NL_MAX_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->max_iterations) ;
    } break ;
    case NL_SYMMETRIC: {
        *params = (NLint)(nlCurrentContext->symmetric) ;
    } break ;
    case NL_USED_ITERATIONS: {
        *params = (NLint)(nlCurrentContext->used_iterations) ;
    } break ;
    case NL_PRECONDITIONER: {
        *params = (NLint)(nlCurrentContext->preconditioner) ;        
    } break ;
    case NL_NNZ: {
        *params = (NLint)(nlSparseMatrixNNZ(&(nlCurrentContext->M)));
    } break;
    default: {
        nlError("nlGetIntegerv","Invalid parameter");
        nl_assert_not_reached ;
    } break ;
    }
}


/* Enable / Disable */

void nlEnable(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        nl_assert(nlCurrentContext->state != NL_STATE_ROW) ;
        nlCurrentContext->normalize_rows = NL_TRUE ;
    } break ;
    case NL_VERBOSE: {
        nlCurrentContext->verbose = NL_TRUE ;
    } break;
    default: {
        nlError("nlEnable","Invalid parameter");        
        nl_assert_not_reached ;
    }
    }
}

void nlDisable(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        nl_assert(nlCurrentContext->state != NL_STATE_ROW) ;
        nlCurrentContext->normalize_rows = NL_FALSE ;
    } break ;
    case NL_VERBOSE: {
        nlCurrentContext->verbose = NL_FALSE ;
    } break;
    default: {
        nlError("nlDisable","Invalid parameter");                
        nl_assert_not_reached ;
    }
    }
}

NLboolean nlIsEnabled(NLenum pname) {
    switch(pname) {
    case NL_NORMALIZE_ROWS: {
        return nlCurrentContext->normalize_rows ;
    } break ;
    case NL_VERBOSE: {
        return nlCurrentContext->verbose ;
    } break ;
    default: {
        nlError("nlIsEnables","Invalid parameter");
        nl_assert_not_reached ;
    }
    }
    return NL_FALSE ;
}


/* NL functions */

void  nlSetFunction(NLenum pname, NLfunc param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        nlCurrentContext->solver_func = nl_cast(NLSolverFunc,param);
        break ;
    case NL_FUNC_MATRIX:
        nlCurrentContext->matrix_vector_prod = nl_cast(NLMatrixFunc,param) ;
        nlCurrentContext->solver = NL_SOLVER_USER ;
        break ;
    case NL_FUNC_PRECONDITIONER:
        nlCurrentContext->precond_vector_prod = nl_cast(NLMatrixFunc,param) ;
        nlCurrentContext->preconditioner = NL_PRECOND_USER ;
        break ;
    case NL_FUNC_PROGRESS:
        nlCurrentContext->progress_func = nl_cast(NLProgressFunc,param) ;
        break ;
    default:
        nlError("nlSetFunction","Invalid parameter");        
        nl_assert_not_reached ;
    }
}

void nlGetFunction(NLenum pname, NLfunc* param) {
    switch(pname) {
    case NL_FUNC_SOLVER:
        *param = nl_cast(NLfunc,nlCurrentContext->solver_func) ;
        break ;
    case NL_FUNC_MATRIX:
        *param = nl_cast(NLfunc,nlCurrentContext->matrix_vector_prod) ;
        break ;
    case NL_FUNC_PRECONDITIONER:
        *param = nl_cast(NLfunc,nlCurrentContext->precond_vector_prod) ;
        break ;
    default:
        nlError("nlGetFunction","Invalid parameter");                
        nl_assert_not_reached ;
    }
}


/* Get/Set Lock/Unlock variables */

void nlSetVariable(NLuint index, NLdouble value) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].value = value ;    
}

NLdouble nlGetVariable(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    return nlCurrentContext->variable[index].value ;
}

void nlLockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].locked = NL_TRUE ;
}

void nlUnlockVariable(NLuint index) {
    nlCheckState(NL_STATE_SYSTEM) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    nlCurrentContext->variable[index].locked = NL_FALSE ;
}

NLboolean nlVariableIsLocked(NLuint index) {
    nl_assert(nlCurrentContext->state != NL_STATE_INITIAL) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    return nlCurrentContext->variable[index].locked  ;
}


/* System construction */

void nlVariablesToVector() {
    NLuint i ;
    nl_assert(nlCurrentContext->alloc_x) ;
    nl_assert(nlCurrentContext->alloc_variable) ;
    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        NLVariable* v = &(nlCurrentContext->variable[i]) ;
        if(!v->locked) {
            nl_assert(v->index < nlCurrentContext->n) ;
            nlCurrentContext->x[v->index] = v->value ;
        }
    }
}

void nlVectorToVariables() {
    NLuint i ;
    nl_assert(nlCurrentContext->alloc_x) ;
    nl_assert(nlCurrentContext->alloc_variable) ;
    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        NLVariable* v = &(nlCurrentContext->variable[i]) ;
        if(!v->locked) {
            nl_assert(v->index < nlCurrentContext->n) ;
            v->value = nlCurrentContext->x[v->index] ;
        }
    }
}


void nlBeginSystem() {
    nlTransition(NL_STATE_INITIAL, NL_STATE_SYSTEM) ;
    nl_assert(nlCurrentContext->nb_variables > 0) ;
    nlCurrentContext->variable = NL_NEW_ARRAY(
        NLVariable, nlCurrentContext->nb_variables
    ) ;
    nlCurrentContext->alloc_variable = NL_TRUE ;
}

void nlEndSystem() {
    nlTransition(NL_STATE_MATRIX_CONSTRUCTED, NL_STATE_SYSTEM_CONSTRUCTED) ;    
}

void nlBeginMatrix() {
    NLuint i ;
    NLuint n = 0 ;
    NLenum storage = NL_MATRIX_STORE_ROWS ;

    
    nlTransition(NL_STATE_SYSTEM, NL_STATE_MATRIX) ;

    for(i=0; i<nlCurrentContext->nb_variables; i++) {
        if(!nlCurrentContext->variable[i].locked) {
            nlCurrentContext->variable[i].index = n ;
            n++ ;
        } else {
            nlCurrentContext->variable[i].index = (NLuint)~0 ;
        }
    }

    nlCurrentContext->n = n ;

    /*
     * If the user trusts OpenNL and has left solver as NL_SOLVER_DEFAULT,
     * then we setup reasonable parameters for him.
     */
    if(nlCurrentContext->solver == NL_SOLVER_DEFAULT) {
        if(nlCurrentContext->least_squares || nlCurrentContext->symmetric) {
            nlCurrentContext->solver = NL_CG;
            if(!nlCurrentContext->preconditioner_defined) {
                nlCurrentContext->preconditioner = NL_PRECOND_JACOBI;
            }
        } else {
            nlCurrentContext->solver = NL_BICGSTAB;
        }
        if(!nlCurrentContext->max_iterations_defined) {
            nlCurrentContext->max_iterations = n*5;
        }
        if(!nlCurrentContext->threshold_defined) {
            nlCurrentContext->threshold = 1e-6;
        }
    }

    
    /* SSOR preconditioner requires rows and columns */
    if(nlCurrentContext->preconditioner == NL_PRECOND_SSOR) {
        storage = (storage | NL_MATRIX_STORE_COLUMNS) ;
    }

    /* a least squares problem results in a symmetric matrix */
    if(
        nlCurrentContext->least_squares  &&
       !nlSolverIsCNC((NLint)(nlCurrentContext->solver))
    ) {
        nlCurrentContext->symmetric = NL_TRUE ;
    }

    if(nlCurrentContext->symmetric) {
        storage = (storage | NL_MATRIX_STORE_SYMMETRIC) ;
    }

    /* SuperLU storage does not support symmetric storage */
    if(
        nlCurrentContext->solver == NL_SUPERLU_EXT       ||
        nlCurrentContext->solver == NL_PERM_SUPERLU_EXT  ||
        nlCurrentContext->solver == NL_SYMMETRIC_SUPERLU_EXT 
    ) {
        storage = (storage & (NLenum)(~NL_SYMMETRIC)) ;
    }

    /*
     * I deactivate symmetric storage, so that the solver
     * can work in parallel mode (see nlSparseMatrix_mult_rows())
     */
    if(nlCurrentContext->preconditioner != NL_PRECOND_SSOR) {
        storage = (storage & (NLenum)(~NL_SYMMETRIC)) ;
    }
    
    nlSparseMatrixConstruct(&nlCurrentContext->M, n, n, storage) ;
    nlCurrentContext->alloc_M = NL_TRUE ;

    nlCurrentContext->x = NL_NEW_ARRAY(NLdouble, n) ;
    nlCurrentContext->alloc_x = NL_TRUE ;
    
    nlCurrentContext->b = NL_NEW_ARRAY(NLdouble, n) ;
    nlCurrentContext->alloc_b = NL_TRUE ;

    nlVariablesToVector() ;

    nlRowColumnConstruct(&nlCurrentContext->af) ;
    nlCurrentContext->alloc_af = NL_TRUE ;
    nlRowColumnConstruct(&nlCurrentContext->al) ;
    nlCurrentContext->alloc_al = NL_TRUE ;
    nlRowColumnConstruct(&nlCurrentContext->xl) ;
    nlCurrentContext->alloc_xl = NL_TRUE ;

    nlCurrentContext->current_row = 0 ;
}

void nlEndMatrix() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_MATRIX_CONSTRUCTED) ;    
    
    nlRowColumnDestroy(&nlCurrentContext->af) ;
    nlCurrentContext->alloc_af = NL_FALSE ;
    nlRowColumnDestroy(&nlCurrentContext->al) ;
    nlCurrentContext->alloc_al = NL_FALSE ;
    nlRowColumnDestroy(&nlCurrentContext->xl) ;
    nlCurrentContext->alloc_al = NL_FALSE ;
    
    if(!nlCurrentContext->least_squares) {
        nl_assert(
            nlCurrentContext->ij_coefficient_called || (
                nlCurrentContext->current_row == 
                nlCurrentContext->n
            )
        ) ;
    }

    nlSparseMatrixComputeDiagInv(&nlCurrentContext->M);
    if(
        nlCurrentContext->preconditioner != NL_PRECOND_SSOR &&
        nlCurrentContext->solver != NL_SUPERLU_EXT       &&
        nlCurrentContext->solver != NL_PERM_SUPERLU_EXT  &&
        nlCurrentContext->solver != NL_SYMMETRIC_SUPERLU_EXT 
       
    ) {
        if(getenv("NL_LOW_MEM") == NULL) {
            nlSparseMatrixCompress(&nlCurrentContext->M);
        }
    }
}

void nlBeginRow() {
    nlTransition(NL_STATE_MATRIX, NL_STATE_ROW) ;
    nlRowColumnZero(&nlCurrentContext->af) ;
    nlRowColumnZero(&nlCurrentContext->al) ;
    nlRowColumnZero(&nlCurrentContext->xl) ;
}

void nlScaleRow(NLdouble s) {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLuint nf            = af->size ;
    NLuint nl            = al->size ;
    NLuint i ;
    for(i=0; i<nf; i++) {
        af->coeff[i].value *= s ;
    }
    for(i=0; i<nl; i++) {
        al->coeff[i].value *= s ;
    }
    nlCurrentContext->right_hand_side *= s ;
}

void nlNormalizeRow(NLdouble weight) {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLuint nf            = af->size ;
    NLuint nl            = al->size ;
    NLuint i ;
    NLdouble norm = 0.0 ;
    for(i=0; i<nf; i++) {
        norm += af->coeff[i].value * af->coeff[i].value ;
    }
    for(i=0; i<nl; i++) {
        norm += al->coeff[i].value * al->coeff[i].value ;
    }
    norm = sqrt(norm) ;
    nlScaleRow(weight / norm) ;
}

void nlEndRow() {
    NLRowColumn*    af = &nlCurrentContext->af ;
    NLRowColumn*    al = &nlCurrentContext->al ;
    NLRowColumn*    xl = &nlCurrentContext->xl ;
    NLSparseMatrix* M  = &nlCurrentContext->M  ;
    NLdouble* b        = nlCurrentContext->b ;
    NLuint nf          = af->size ;
    NLuint nl          = al->size ;
    NLuint current_row = nlCurrentContext->current_row ;
    NLuint i ;
    NLuint j ;
    NLdouble S ;
    nlTransition(NL_STATE_ROW, NL_STATE_MATRIX) ;

    if(nlCurrentContext->normalize_rows) {
        nlNormalizeRow(nlCurrentContext->row_scaling) ;
    } else {
        nlScaleRow(nlCurrentContext->row_scaling) ;
    }
    /*
     * if least_squares : we want to solve
     * A'A x = A'b
     */
    if(nlCurrentContext->least_squares) {
        for(i=0; i<nf; i++) {
            for(j=0; j<nf; j++) {
                nlSparseMatrixAdd(
                    M, af->coeff[i].index, af->coeff[j].index,
                    af->coeff[i].value * af->coeff[j].value
                ) ;
            }
        }
        S = -nlCurrentContext->right_hand_side ;
        for(j=0; j<nl; j++) {
            S += al->coeff[j].value * xl->coeff[j].value ;
        }
        for(i=0; i<nf; i++) {
            b[ af->coeff[i].index ] -= af->coeff[i].value * S ;
        }
    } else {
        for(i=0; i<nf; i++) {
            nlSparseMatrixAdd(
                M, current_row, af->coeff[i].index, af->coeff[i].value
            ) ;
        }
        b[current_row] = nlCurrentContext->right_hand_side ;
        for(i=0; i<nl; i++) {
            b[current_row] -= al->coeff[i].value * xl->coeff[i].value ;
        }
    }
    nlCurrentContext->current_row++ ;
    nlCurrentContext->right_hand_side     = 0.0 ;
    nlCurrentContext->right_hand_side_set = NL_FALSE ;
    nlCurrentContext->row_scaling         = 1.0 ;
}

void nlCoefficient(NLuint index, NLdouble value) {
    NLVariable* v = NULL ;
    nlCheckState(NL_STATE_ROW) ;
    nl_debug_range_assert(index, 0, nlCurrentContext->nb_variables - 1) ;
    v = &(nlCurrentContext->variable[index]) ;
    if(v->locked) {
        nlRowColumnAppend(&(nlCurrentContext->al), 0, value) ;
        nlRowColumnAppend(&(nlCurrentContext->xl), 0, v->value) ;
    } else {
        nlRowColumnAppend(&(nlCurrentContext->af), v->index, value) ;
    }
}

void nlAddIJCoefficient(NLuint i, NLuint j, NLdouble value) {
    NLSparseMatrix* M  = &nlCurrentContext->M;    
    nlCheckState(NL_STATE_MATRIX);
    nl_debug_range_assert(i, 0, nlCurrentContext->nb_variables - 1);
    nl_debug_range_assert(j, 0, nlCurrentContext->nb_variables - 1);
#ifdef NL_DEBUG
    for(NLuint i=0; i<nlCurrentContext->nb_variables; ++i) {
        nl_debug_assert(!nlCurrentContext->variable[i].is_locked);
    }
#endif    
    nlSparseMatrixAdd(M, i, j, value);
    nlCurrentContext->ij_coefficient_called = NL_TRUE;
}

void nlAddIRightHandSide(NLuint i, NLdouble value) {
    nlCheckState(NL_STATE_MATRIX);
    nl_debug_range_assert(i, 0, nlCurrentContext->nb_variables - 1);
#ifdef NL_DEBUG
    for(NLuint i=0; i<nlCurrentContext->nb_variables; ++i) {
        nl_debug_assert(!nlCurrentContext->variable[i].is_locked);
    }
#endif
    nlCurrentContext->b[i] += value;
    nlCurrentContext->ij_coefficient_called = NL_TRUE;
}


void nlRightHandSide(NLdouble value) {
    nl_assert(!nlCurrentContext->right_hand_side_set);
    nlCurrentContext->right_hand_side = value;
    nlCurrentContext->right_hand_side_set = NL_TRUE;
}

void nlRowScaling(NLdouble value) {
    nlCheckState(NL_STATE_MATRIX);
    nlCurrentContext->row_scaling = value;
}

void nlBegin(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlBeginSystem() ;
    } break ;
    case NL_MATRIX: {
        nlBeginMatrix() ;
    } break ;
    case NL_ROW: {
        nlBeginRow() ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}

void nlEnd(NLenum prim) {
    switch(prim) {
    case NL_SYSTEM: {
        nlEndSystem() ;
    } break ;
    case NL_MATRIX: {
        nlEndMatrix() ;
    } break ;
    case NL_ROW: {
        nlEndRow() ;
    } break ;
    default: {
        nl_assert_not_reached ;
    }
    }
}


/* nlSolve() driver routine */

NLboolean nlSolve() {
    NLboolean result ;
    NLdouble start_time = nlCurrentTime() ; 
    nlCheckState(NL_STATE_SYSTEM_CONSTRUCTED) ;
    nlCurrentContext->elapsed_time = 0.0 ;
    nlCurrentContext->flops = 0 ;    
    result =  nlCurrentContext->solver_func() ;
    nlVectorToVariables() ;
    nlCurrentContext->elapsed_time = nlCurrentTime() - start_time ;
    nlTransition(NL_STATE_SYSTEM_CONSTRUCTED, NL_STATE_SOLVED) ;
    return result ;
}

