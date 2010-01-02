/*!
 * \file
 * \brief Definitions of functions for solving linear equation systems
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef LS_SOLVE_H
#define LS_SOLVE_H

#include <itpp/base/mat.h>


namespace itpp
{


/*! \addtogroup linearequations
 */
//!@{

/*! \brief Solve linear equation system by LU factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a \f$n \times n\f$ matrix.
Uses the LAPACK routine DGESV.
*/
bool ls_solve(const mat &A, const vec &b, vec &x);

/*! \brief Solve linear equation system by LU factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a \f$n \times n\f$ matrix.
Uses the LAPACK routine DGESV.
*/
vec ls_solve(const mat &A, const vec &b);

/*! \brief Solve multiple linear equations by LU factorisation.

Solves the linear system \f$AX=B\f$. Here \f$A\f$ is a nonsingular \f$n \times n\f$ matrix.
Uses the LAPACK routine DGESV.
*/
bool ls_solve(const mat &A, const mat &B, mat &X);

/*! \brief Solve multiple linear equations by LU factorisation.

Solves the linear system \f$AX=B\f$. Here \f$A\f$ is a nonsingular \f$n \times n\f$ matrix.
Uses the LAPACK routine DGESV.
*/
mat ls_solve(const mat &A, const mat &B);


/*! \brief Solve linear equation system by LU factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a \f$n \times n\f$ matrix.
Uses the LAPACK routine ZGESV.
*/
bool ls_solve(const cmat &A, const cvec &b, cvec &x);

/*! \brief Solve linear equation system by LU factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a \f$n \times n\f$ matrix.
Uses the LAPACK routine ZGESV.
*/
cvec ls_solve(const cmat &A, const cvec &b);

/*! \brief Solve multiple linear equations by LU factorisation.

Solves the linear system \f$AX=B\f$. Here \f$A\f$ is a nonsingular \f$n \times n\f$ matrix.
Uses the LAPACK routine ZGESV.
*/
bool ls_solve(const cmat &A, const cmat &B, cmat &X);

/*! \brief Solve multiple linear equations by LU factorisation.

Solves the linear system \f$AX=B\f$. Here \f$A\f$ is a nonsingular \f$n \times n\f$ matrix.
Uses the LAPACK routine ZGESV.
*/
cmat ls_solve(const cmat &A, const cmat &B);


/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a symmetric positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine DPOSV.
*/
bool ls_solve_chol(const mat &A, const vec &b, vec &x);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a symmetric positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine DPOSV.
*/
vec ls_solve_chol(const mat &A, const vec &b);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$AX=B\f$, where \f$A\f$ is a symmetric positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine DPOSV.
*/
bool ls_solve_chol(const mat &A, const mat &B, mat &X);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$AX=B\f$, where \f$A\f$ is a symmetric positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine DPOSV.
*/
mat ls_solve_chol(const mat &A, const mat &B);


/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a Hermitian positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine ZPOSV.
*/
bool ls_solve_chol(const cmat &A, const cvec &b, cvec &x);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$Ax=b\f$, where \f$A\f$ is a Hermitian positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine ZPOSV.
*/
cvec ls_solve_chol(const cmat &A, const cvec &b);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$AX=B\f$, where \f$A\f$ is a Hermitian positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine ZPOSV.
*/
bool ls_solve_chol(const cmat &A, const cmat &B, cmat &X);

/*! \brief Solve linear equation system by Cholesky factorisation.

Solves the linear system \f$AX=B\f$, where \f$A\f$ is a Hermitian positive definite \f$n \times n\f$ matrix.
Uses the LAPACK routine ZPOSV.
*/
cmat ls_solve_chol(const cmat &A, const cmat &B);



/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and is built upon the LAPACK routine DGELS.
*/
bool ls_solve_od(const mat &A, const vec &b, vec &x);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
vec ls_solve_od(const mat &A, const vec &b);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
bool ls_solve_od(const mat &A, const mat &B, mat &X);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
mat ls_solve_od(const mat &A, const mat &B);


/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and is built upon the LAPACK routine ZGELS.
*/
bool ls_solve_od(const cmat &A, const cvec &b, cvec &x);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
cvec ls_solve_od(const cmat &A, const cvec &b);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
bool ls_solve_od(const cmat &A, const cmat &B, cmat &X);

/*! \brief Solves overdetermined linear equation systems.

Solves the overdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \geq n\f$.
Uses QR-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
cmat ls_solve_od(const cmat &A, const cmat &B);



/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and is built upon the LAPACK routine DGELS.
*/
bool ls_solve_ud(const mat &A, const vec &b, vec &x);

/*! \brief Solves overdetermined linear equation systems.

Solves the underdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
vec ls_solve_ud(const mat &A, const vec &b);

/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
bool ls_solve_ud(const mat &A, const mat &B, mat &X);

/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine DGELS.
*/
mat ls_solve_ud(const mat &A, const mat &B);


/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and is built upon the LAPACK routine ZGELS.
*/
bool ls_solve_ud(const cmat &A, const cvec &b, cvec &x);

/*! \brief Solves overdetermined linear equation systems.

Solves the underdetermined linear system \f$Ax=b\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
cvec ls_solve_ud(const cmat &A, const cvec &b);

/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
bool ls_solve_ud(const cmat &A, const cmat &B, cmat &X);

/*! \brief Solves underdetermined linear equation systems.

Solves the underdetermined linear system \f$AX=B\f$, where \f$A\f$ is a \f$m \times n\f$ matrix and \f$m \leq n\f$.
Uses LQ-factorization and assumes that \f$A\f$ is full rank. Based on the LAPACK routine ZGELS.
*/
cmat ls_solve_ud(const cmat &A, const cmat &B);


/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,b,x), ls_solve_od(A,b,x) or ls_solve_ud(A,b,x)
*/
bool backslash(const mat &A, const vec &b, vec &x);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,b), ls_solve_od(A,b) or ls_solve_ud(A,b)
*/
vec backslash(const mat &A, const vec &b);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,B,X), ls_solve_od(A,B,X), or ls_solve_ud(A,B,X).
*/
bool backslash(const mat &A, const mat &B, mat &X);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,B), ls_solve_od(A,B), or ls_solve_ud(A,B).
*/
mat backslash(const mat &A, const mat &B);


/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,b,x), ls_solve_od(A,b,x) or ls_solve_ud(A,b,x)
*/
bool backslash(const cmat &A, const cvec &b, cvec &x);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,b), ls_solve_od(A,b) or ls_solve_ud(A,b)
*/
cvec backslash(const cmat &A, const cvec &b);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,B,X), ls_solve_od(A,B,X), or ls_solve_ud(A,B,X).
*/
bool backslash(const cmat &A, const cmat &B, cmat &X);

/*! \brief A general linear equation system solver.

Tries to emulate the backslash operator in Matlab by calling
ls_solve(A,B), ls_solve_od(A,B), or ls_solve_ud(A,B).
*/
cmat backslash(const cmat &A, const cmat &B);



/*! \brief Forward substitution of square matrix.

Solves Lx=b, where L is a lower triangular n by n matrix.
Assumes that L is nonsingular. Requires n^2 flops.
Uses Alg. 3.1.1 in Golub & van Loan "Matrix computations", 3rd ed., p. 89.
*/
vec forward_substitution(const mat &L, const vec &b);

/*! \brief Forward substitution of square matrix.

Solves Lx=b, where L is a lower triangular n by n matrix.
Assumes that L is nonsingular. Requires n^2 flops.
Uses Alg. 3.1.1 in Golub & van Loan "Matrix computations", 3rd ed., p. 89.
*/
void forward_substitution(const mat &L, const vec &b, vec &x);

/*! \brief Forward substitution of band matrices.

Solves Lx=b, where L is a lower triangular n by n band-matrix with lower
bandwidth p.
Assumes that L is nonsingular. Requires about 2np flops (if n >> p).
Uses Alg. 4.3.2 in Golub & van Loan "Matrix computations", 3rd ed., p. 153.
*/
vec forward_substitution(const mat &L, int p, const vec &b);

/*! \brief Forward substitution of band matrices.

Solves Lx=b, where L is a lower triangular n by n band-matrix with
lower bandwidth p.
Assumes that L is nonsingular. Requires about 2np flops (if n >> p).
Uses Alg. 4.3.2 in Golub & van Loan "Matrix computations", 3rd ed., p. 153.
*/
void forward_substitution(const mat &L, int p, const vec &b, vec &x);

/*! \brief Backward substitution of square matrix.

Solves Ux=b, where U is a upper triangular n by n matrix.
Assumes that U is nonsingular. Requires n^2 flops.
Uses Alg. 3.1.2 in Golub & van Loan "Matrix computations", 3rd ed., p. 89.
*/
vec backward_substitution(const mat &U, const vec &b);

/*! \brief Backward substitution of square matrix.

Solves Ux=b, where U is a upper triangular n by n matrix.
Assumes that U is nonsingular. Requires n^2 flops.
Uses Alg. 3.1.2 in Golub & van Loan "Matrix computations", 3rd ed., p. 89.
*/
void backward_substitution(const mat &U, const vec &b, vec &x);

/*! \brief Backward substitution of band matrix.

Solves Ux=b, where U is a upper triangular n by n matrix band-matrix with
upper bandwidth q.
Assumes that U is nonsingular. Requires about 2nq flops (if n >> q).
Uses Alg. 4.3.3 in Golub & van Loan "Matrix computations", 3rd ed., p. 153.
*/
vec backward_substitution(const mat &U, int q, const vec &b);

/*! \brief Backward substitution of band matrix.

Solves Ux=b, where U is a upper triangular n by n matrix band-matrix with
upper bandwidth q.
Assumes that U is nonsingular. Requires about 2nq flops (if n >> q).
Uses Alg. 4.3.3 in Golub & van Loan "Matrix computations", 3rd ed., p. 153.
*/
void backward_substitution(const mat &U, int q, const vec &b, vec &x);

//!@}

} //namespace itpp

#endif // #ifndef LS_SOLVE_H



