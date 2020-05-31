#lang scribble/manual
@;raco scribble +m --redirect-main http://docs.racket-lang.org flmatrix.scrbl && open flmatrix.html
@(require scribble/example)
@(define reference.scrbl '(lib "scribblings/reference/reference.scrbl"))
@(define math.scrbl      '(lib "math/scribblings/math.scrbl"))
@(require (for-label racket/base ffi/vector ffi/unsafe))

@title[#:tag "flmatrix"]{Flmatrix: Floating Point Matrices}

@defmodule[flmatrix]

This manual documents the matrix library @racketmodname[flmatrix].

@author[@author+email["Jens Axel Søgaard" "jensaxel@soegaard.net"]]

@local-table-of-contents[#:style 'immediate-only]

@section{Introduction}

A matrix is a rectangular arrangements of numbers in rows and columns.
This library provides functions to construct and compute with
matrices whose elements are IEEE double precision floating point numbers.
These numbers are referred to as @tech[#:doc reference.scrbl]{flonums}
in the Racket manual, but the most common name for these numbers are 
simply @emph{doubles}.

Restricting the scope of the library to matrices with floating numbers
allow the implementation to use routines implemented in Fortran and C.
The low-level routines consists of calls to functions in BLAS and LAPACK.
BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra PACKage)
are industry standard libraries and are available on all major platforms.

If you are in need of matrix algebra over more general numbers then
look at the functional matrix library in @secref["matrices" #:doc math.scrbl].

This library can be used in a functional manner, but imperative operations
are available. There are at least two reasons to use the imperative approach:
1) text books often describe matrix algorithms using imperative operations,
and, 2) imperative operations can reduce the amount of memory needed
during a computation.

The available operations can be divided into rough categories:
@itemlist[
  @item{Level 1:  High   level - do what I mean}
  @item{Level 2:  Medium level - do what I mean this way}
  @item{Level 3:  Low    level - do it using this underlying C-function}]

To use the library one can get by with level 1 operations, but if you understand
the underlying representation, you can improve your algorithms using 
level 2 operations. For those that really want to squeeze out the last bit of
performance we have made level 3 operations available as well.


@section{Quick Tutorial}

@(define quick-eval (let ([e (make-base-eval)]) (e '(require flmatrix/flmatrix)) e))

This section shows how to do simple matrix computations.
The beginning part of the tutorial describes working with matrices simply as arrays of numbers.
The end shows how to do linear algebra.

@subsection{Basic Properties}

An @racket[flmatrix] consists conceptually of a two-dimensional
array of floating point numbers. An mxn (m by n) matrix is
divided into m rows and n columns. 

The basic properties of an @racket[flmatrix] can be examined using these functions:

@bold[@racket[(shape A)]]
  return a list of with the number of rows and columns

@bold[@racket[(size A)]]
  the number of elements in the matrix

@bold[@racket[(nrows A)]]
  the number of rows

@bold[@racket[(ncols A)]]
  the number of columns

@examples[#:eval quick-eval
          (define A (flmatrix: [[1 2 3]
                                [4 5 5]]))
          (shape A)   ; dimensions 
          (size  A)
          (nrows A)
          (ncols A)]

@subsection{Basic Indexing}

Since a matrix is divided into rows and columns we can refer to an
element in the matrix by row and column numbers. The element on the
i'th row and j'th column is referred to as the element with index (i,j).

The indices are zero-based so a matrix with m rows and n columns
has row-indices 0, 1, ..., m-1 and column-indices 0, 1, ... n-1.

@bold[@racket[(ref A i j)]]
  the element in A with index (i,j)

@bold[@racket[(row A i)]]
  the i'th row of A

@bold[@racket[(col A j)]]
  the j'th row of A

Notice that row and column vectors are simply matrices with
a single row and a single column respectively

@examples[#:eval quick-eval
          (define A (flmatrix: [[1 2 3]
                                [4 5 5]]))
          (ref A 0 1)
          (row A 0)
          (col A 1)]


@subsection{Matrix Creation}

There are several ways of creating matrices.

Use @racket[matrix] to create an @racket[flmatrix] from existing Racket data.
It can convert vector-of-vector-of and list-of-lists representation of matrices
into the @racket[flmatrix] representation. A vector of numbers or a list
of numbers will be converted into a column vector (a matrix with only one column).

Any non-floating point numbers will be converted to floating point.
The function @racket[matrix] also accepts @racket[f64vectors] as input.

@bold[@racket[(matrix obj)]]
  create a matrix with values from @racket[obj]

@examples[#:eval quick-eval
          (matrix '[[1/2 1/3] [4 5]])
          (matrix #[#[1 2 3] #[4 5 6]])
          (matrix (list 1 2 3))
          (matrix (vector 1 2 3))
          (matrix (f64vector 1 2 3))]

After conversion the created @racket[flmatrix] will contain a pointer to
piece newly allocated memory containing the floating point numbers.
If you happen to work with data in the form of @racket[f64vector]s, then
you can avoid the allocation, if you use @racket[matrix!] instead.
If the same @racket[f64vector] is used to create two matrices with @racket[matrix!]
they will share the same backing array - so setting an element one matrix
will affect the other.

@bold[@racket[(matrix! obj)]]
  create a matrix with values from @racket[obj] avoid allocation of
  backing array if possible

@examples[#:eval quick-eval
          (define v (f64vector 1 2 3))
          (define A (matrix! v))
          (define B (matrix! v))
          (list A B)
          (mset! A 0 0 42)
          (list A B)]
For comparision the same example with @racket[matrix] instead:
@examples[#:eval quick-eval
          (define v (f64vector 1 2 3))
          (define A (matrix v))
          (define B (matrix v))
          (list A B)
          (mset! A 0 0 42)
          (list A B)]

In order to create a matrix of specific size with all zeros or all ones,
use the functions @racket[zeros] and @racket[ones].


@bold[@racket[(zeros n)]]
  create a square nxn matrix with all zeros

@bold[@racket[(zeros m n)]]
  create a mxn matrix with all zeros

@bold[@racket[(ones n)]]
  create a square nxn matrix with all ones

@bold[@racket[(ones m n)]]
  create a mxn matrix with all ones

@examples[#:eval quick-eval
          (zeros 2)
          (zeros 2 3)
          (ones 2)
          (ones 2 3)]

To create ranges of values use @racket[arange] or @racket[colarange] which both work like
@racket[(matrix (range start stop step))], but avoids build an intermediary list.
The functions @racket[arange] and @racket[colarange] produce row and column vectors respectively.
The vector created has length @racket[(ceiling (/ (- stop start) step))].

@bold[@racket[(arange start stop step)]]
  create a row vector with values from start to stop (exclusively),
  here step is the gap between values

@bold[@racket[(arange start stop)]]
  like @racket[(arange start stop 1.0)]

@bold[@racket[(arange stop)]]
  like @racket[(arange 0.0 stop 1.0)]

@bold[@racket[(colarange start stop step)]] @linebreak[]
@bold[@racket[(colarange start stop)]]      @linebreak[]
@bold[@racket[(colarange start)]]           @linebreak[]
  like @racket[arange] but produces a column vector.


@examples[#:eval quick-eval
          (arange 5 10 2)
          (arange 5 10)
          (arange 5)]
@examples[#:eval quick-eval
          (colarange 5 10 2)
          (colarange 5 10)
          (colarange 5)]

Sometimes it is possible to keep the elements of matrix, but change its shape.

@bold[@racket[(reshape A m n)]]  @linebreak[]
  return a matrix with shape mxn using the elements of A,

@bold[@racket[(reshape! A m n)]] @linebreak[]
  return a matrix with shape mxn using the elements of A, share the backing area with A

@examples[#:eval quick-eval
          (arange 9)
          (reshape (arange 9) 3 3)
          (transpose (reshape (arange 9) 3 3))]

As an alternative to @racket[arange] consider using @racket[linspace], which
allow you to provide an exact endpoint.

@bold[@racket[(linspace start stop num)]] @linebreak[]
  return a column vector with @racket[num] numbers evenly spaced from
  @racket[start] to @racket[stop]

@bold[@racket[(linspace start stop num #f)]] @linebreak[]
   like @racket[(linspace start stop num)] but omit the last number

@examples[#:eval quick-eval
          (linspace 2 4 6)
          (linspace 2 4 6 #f)]


@subsection{Elementwise Operations}

Elementwise operations (also called @emph{pointwise} operations) work on each element.
The operations are named with a beginning point. Beside the elementwise versions
of the standard arithmetic operations, the standard numerical functions also
has an elementwise counterpart. Binary operators work both on matrices (of the same side)
and on a number and matrix.

@bold[@racket[(.+ A B)]]                            @linebreak[]
@bold[@racket[(.- A B)]] and @bold[@racket[(.- A)]] @linebreak[]
@bold[@racket[(.* A B)]] @linebreak[]
@bold[@racket[(./ A B)]] and @bold[@racket[(./ A)]] @linebreak[]
  Elementwise version of the arithmetical operations. 
  The operations returns the result as a new matrix.

Note that @racket[.*] is elementwise multiplication. Use @racket[times]
to multiply two matrices in the linear algebra sense.

@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.- A)
          (./ A)
          (.+ A B)
          (.- A B)
          (.* A B)
          (./ A B)]
One of the arguments can be a number:
@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.+ A 1)
          (.- A 2)
          (.* 3 B)
          (./ A 4)]

The elementwise versions of the standard numerical functions are:

@bold[@racket[(.sin A)]]
@bold[@racket[(.cos A)]]
@bold[@racket[(.tan A)]]
@bold[@racket[(.exp A)]]
@bold[@racket[(.log A)]]
@bold[@racket[(.sqr A)]]
@bold[@racket[(.sqrt A)]]
@bold[@racket[(.expt A B)]]

@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.sqr A)
          (.expt A 2)
          (.expt 2 A)]


The elementwise operations above all, allocate a new matrix.
If instead you want to modify the elements of an existing matrix,
the following functions are for you.

@bold[@racket[(.-! A)]]
@bold[@racket[(./! A)]]
@bold[@racket[(.sin! A)]]
@bold[@racket[(.cos! A)]]
@bold[@racket[(.tan! A)]]
@bold[@racket[(.exp! A)]]
@bold[@racket[(.log! A)]]
@bold[@racket[(.sqr! A)]]
@bold[@racket[(.sqrt! A)]]
@bold[@racket[(.expt! A B)]]
@bold[@racket[(.+! A B)]]      
@bold[@racket[(.-! A B)]]
@bold[@racket[(.*! A B)]]
@bold[@racket[(./! A B)]] @linebreak[]


For binary operations, the result is stored in the first argument.
@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (.-! A)
          (.-! A)
          (.expt! A 2)
          (.expt! A 2)]

Also, if you want to store the result of an elementwise in another 
matrix C, you can do as follows for the unary operations:


@bold[@racket[(.sin!  A C)]]
@bold[@racket[(.cos!  A C)]]
@bold[@racket[(.tan!  A C)]]
@bold[@racket[(.exp!  A C)]]
@bold[@racket[(.log!  A C)]]
@bold[@racket[(.sqr!  A C)]]
@bold[@racket[(.sqrt! A C)]]

And for the binary operations:

@bold[@racket[(.+! A B C)]]      
@bold[@racket[(.-! A B C)]]
@bold[@racket[(.*! A B C)]]
@bold[@racket[(./! A B C)]] @linebreak[]
@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.sqr! B A)
          A]

Finally, for @racket[.-!] and @racket[./!] which are both unary and binary operations
at once, use @racket[#f] as @racket[B] to get the unary version.

@examples[#:eval quick-eval
          (define A (matrix '((0 1) (2 3))))
          (define B (matrix '((4 5) (6 7))))
          (.-! B #f A)
          A]


@subsection{Indexing, Submatrices and Iterating}

From the section on Basic Indexing we know that the element on row i in column j,
has index (i,j) and can be extracted with the function @racket[ref].

The i'th row and the j'th column can be extraced with @racket[row] and @racket[col]
respectively.

To get a submatrix use @racket[sub] and @racket[sub!].

@bold[@racket[(sub A i j m n)]] @linebreak[]
  Make a copy of the submatrix of A with upper left corner in (i,j) and with size mxn.

@bold[@racket[(sub! A i j m n)]] @linebreak[]
  Same as @racket[sub], but the elements are copied - the underlying 
  array of flonums are shared.


@examples[#:eval quick-eval
          (define A (transpose (reshape (arange 25) 5 5)))
          A
          (sub A 0 0 3 2)
          (sub A 1 1 2 2)]

The function @racket[sub!] can be used to mutate part of a larger submatrix.

Let's say we have a matrix, in which we want to zero out all elements except
those on the edges. We can use @racket[sub!] to get a submatrix of the inner part,
then use @racket[zeros!] to clear the elements.

@examples[#:eval quick-eval
          (define A (transpose (reshape (arange 10 35) 5 5)))
          A
          (define B (sub! A 1 1 3 3))
          B
          (zeros! B)
          A]

To iterate over a row or a column use @racket[in-flrow] and @racket[in-flcolumn].

@examples[#:eval quick-eval
          (define A (matrix '((11 22) (33 44))))
          (for/list ([   x  (in-flrow A 0)]) x)
          (for/list ([(i x) (in-flrow A 0)]) (list x i))
          (for/list ([   x  (in-flcolumn A 0)]) x)
          (for/list ([(i x) (in-flcolumn A 0)]) (list x i))]

@subsection{Basic Linear Algebra}

The basic linear algebra are @racket[plus], @racket[minus] and @racket[times],
which compute the sum, difference and product of a series of matrices.

@bold[@racket[(plus  A ...)]]
@bold[@racket[(minus A ...)]]
@bold[@racket[(times A ...)]] @linebreak[]
Computes the sum, difference and product of a series of matrices and/or numbers.

@examples[#:eval quick-eval
          (define A (matrix '((2 0) (0 2))))
          (define B (matrix '((1 2) (3 4))))
          (define C (column 4 5))
          (plus A B)
          (plus A 10)
          (plus A 10 B)
          (minus A)
          (minus A B)
          (times A B)
          (times A 2 B)
          (times A C)]

As usual, there are variants that mutate the first given matrix
instead of allocating a new backing array of flonums.

@bold[@racket[(plus!  A B ...)]]
@bold[@racket[(minus! A B ...)]]  @linebreak[]
Like @racket[plus] and @racket[minus]  but stores
the result in @racket[A], which must be a matrix.

@examples[#:eval quick-eval
          (define A (matrix '((2 0) (0 2))))
          (define B (matrix '((0 2) (2 0))))
          (plus! A B)
          A]

@bold[@racket[(power A n)]] @linebreak[]
Computes the n'th power of a matrix A, where n is a natural number.

@examples[#:eval quick-eval
          (define A (matrix '((1 1) (0 1))))
          (list (power A 0) (power A 1) (power A 2) (power A 3))]

@subsection{Matrix and Vector Products}

The inner product (also known as the dot product) of two column vectors
can be computed by @racket[dot].

@bold[@racket[(dot v w)]] @linebreak[]
Computes the inner product of two column vectors (i.e. matrices with only one column).

@examples[#:eval quick-eval
          (define v (column -1 1))
          (define w (matrix '((2) (2))))
          (dot v w)]

The outer product of a column vector A with m rows and an row B with n columns
is an mxn matrix O with elements o_ij = a_i *b_j.

@bold[@racket[(outer A B)]] @linebreak[]
Computes the outer product of the first column of A and the first row of B.

@examples[#:eval quick-eval
          (define A (column 2 3))
          (define B (transpose (column 5 7)))
          (outer A B)]


The Kronecker product between two matrices @racket[A] and @racket[B] replaces
each element @racket[a] of @racket[A] with a copy of @racket[B] scaled with @racket[A].
The Kronecker product is a generalization of the outer product.

@bold[@racket[(kron A B)]] @linebreak[]
Computes the Kronecker product of the matrices @racket[A] and @racket[B].

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((1 1) (1 1))))
          (kron A B)]


@subsection{Matrix Decompositions}

@bold[@racket[(cholesky A)]] @bold[@racket[(qr A)]] @bold[@racket[(svd A)]] @linebreak[]
Computes the Cholesky, QR and SVD decompositions respectively.


The Singular Value Decomposition (SVD) returns three matrices: 
a unitary matrix U, a column vector of singular values S and 
a unitary matrix VT (V transposed). The function @racket[diag] constructs
a diagonal matrix from the singular values.

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define-values (U S VT) (svd A))
          (define Σ (diag S))
          (list U Σ VT S)
          (times U Σ VT)]

The QR Decomposition of A consists of two matrices: an orthogonal matrix Q
and an upper triangular matrix R such that A=QR.

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define-values (Q R) (qr A))
          (list Q R)
          (times Q R)]

If the matrix A is symmetric and positive-definite, then
the Cholesky decomposition can be computed. 
It comes in two forms

A = L L^T    or    A = U^T U,

where L and U are lower and upper triangular matrices.

Note: @racket[cholesky] does not check that the input matrix @racket[A]
is symmetric and positive definite.

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (2 4))))
          (define L (cholesky A))
          (list L (transpose L))
          (times L (transpose L))
          (define U (cholesky A 'upper))
          (list (transpose U) U)          
          (times (transpose U) U)]


@subsection{Matrix Eigenvalues and Eigenvectors}

Eigenvalues and eigenvectors of a square matrix can be computed with @racket[eig] 
or, if only the eigenvalues are needed, with @racket[eigvals].
Note that even if all elements of a matrix are real, the eigenvalues in some
cases are complex. Therefore the eigenvalues are returned as a standard
Racket vector.


@bold[@racket[(eig A)]] @linebreak[]
Compute eigenvalues and right eigenvectors.

@bold[@racket[(eigvals A)]] @linebreak[]
Compute eigenvalues.

@examples[#:eval quick-eval
          (eig     (diag '(1 2)))
          (eigvals (diag '(1 2)))
          (eig     (matrix '((1 -1) (1 1))))]
          


@subsection{Norms and Invariants}

The standard 2-norm can be computed by @racket[norm].
For a column vector the norm is sometimes referred to as the length.

@bold[@racket[(norm A)]] @linebreak[]
Compute the square root om the sum of the square of all elements.

@examples[#:eval quick-eval
          (norm    (matrix '((1 1))))
          (norm    (matrix '((1 -1) (-1 1))))]


@bold[@racket[(det A)]] @linebreak[]
Computes the determinant of a square matrix @racket[A].

@examples[#:eval quick-eval
          (det  (matrix '((1 2) (0 4))))
          (det  (matrix '((1 1) (2 2))))]

@bold[@racket[(trace A)]] @linebreak[]
Computes the trace, the sum along a diagonal, of a matrix.

@examples[#:eval quick-eval
          (trace (matrix '((1 2) (0 4))))]


@bold[@racket[(rank A)]] @linebreak[]
Computes the rank of a square matrix.
The rank is the dimension of the column space,
which is equal to the dimension of the row space,
which is equal to the number  of non-zero singular values
in an SVD decomposition.

@examples[#:eval quick-eval
          (rank  (matrix '((1 2) (0 4))))
          (rank  (matrix '((1 1) (2 2))))]


@subsection{Solving Equations and Inverting Matrices}

Solving linear equations are more or less the raison d'etre for matrices.
The main workhorse is @racket[mldivide], which can solve for X
in the equation:

    @racket[AX = B],

where @racket[A] is a an mxm matrix, and both X and B are mxn.

Note that @racket[A] needs to be of full rank for the equation
to have a solution. The solver doesn't check that the input
matrix has full rank, it just runs it computation as usual.
To check that the output from @racket[solve] is indeed a solution,
you can evaluate @racket[(times A X)] and compare with @racket[B].
The name @racket[mldivide] is short for "Matrix Left divide".
Although @racket[mldivide] doesn't find @racket[X] by 
multiplying @racket[B] with @racket[A^-1] on the left, 
it is fitting analogy.

@bold[@racket[(mldivide A B)]] @linebreak[]
Solve the equation @racket[AX = B] using LU-decomposition with
partial pivoting. The matrix A must be square and of full rank, the number 
of rows in @racket[A] must be the same as the number columns in @racket[B].

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((1) (0))))
          (define X (mldivide A B))
          (list X (times A X))]

@bold[@racket[(mrdivide B A)]] @linebreak[]
Solve the equation @racket[XA = B].
The name @racket[mrdivide] is short for "Matrix Right divide".

@examples[#:eval quick-eval
          (define A (matrix '((1 2) (3 4))))
          (define B (matrix '((2 4) (6 8))))
          (define X (mrdivide B A))
          (list X (times X A))]


@bold[@racket[(inv A)]] @linebreak[]
Find the multiplicative inverse of a square matrix @racket[A].

@examples[#:eval quick-eval
          (define A    (matrix '((1 2) (3 4))))
          (define Ainv (inv A))
          (list Ainv (times A Ainv))]
An inverse of @racket[A] can be used to solve @racket[AX=B], but
using @racket[mldivide] directly is normally better. However, let's
try to solve the equation from the previous example.

@examples[#:eval quick-eval
          (define B    (matrix '((1) (0))))
          (define X    (times Ainv B))
          (list X (times A X))]


@bold[@racket[(pinv A)]] @linebreak[]
Find the Moore-Penrose pseudo-inverse of the matrix @racket[A].
The matrix @racket[A] does not need to be square.
The pseudo inverse of an mxn matrix is of size 


@examples[#:eval quick-eval
          (define A  (matrix '((1 2) (3 4))))
          (define A+ (pinv A))
          (list A+ (times A+ A A+) (times A A+ A))]
@examples[#:eval quick-eval
          (define B  (matrix '((1 2 3) (4 5 6))))
          (define B+ (pinv B))
          (list B+ (times B+ B B+) (times B B+ B))]


@subsection{Least Square Problems}

Let @racket[A] be an mxn matrix and let @racket[b] be an nx1 column vector.
The equation @racket[Ax=b] may (depending on #racket[A]) may not have
an unique solution - or a solution at all.

As an alternative, one can look for the vector @racket[x] that minimizes:
    @racket[norm(Ax-b)],
where #@racket[norm] is the Euclidean 2-norm.

The function @racket[lstsq] return the minimum norm solution @racket[x]
of the above the problem.

If @racket[lstsq] is given an nxk matrix @racket[B], then the
problem will be solved for each column @racket[b] of @racket[B].



@bold[@racket[(lstsq A B)]] @linebreak[]
Find minimum norm solution to the least squares problem: @racket["minimize |Ax-b|"] ,
for each column b of a larger matrix B.


As an example, let's look at estimating @racket[b0] and @racket[b1] in the model:
    @racket[y=b0*x+b1]
given a data set consisting of corresponding x- and y-values. The calculation reveals
that the relation between @racket[x] and @racket[y] is @racket[y=2x+1].

The matrix @racket[X] is called the @emph{design matrix} of the problem.
See @hyperlink["https://en.wikipedia.org/wiki/Design_matrix"]{Design Matrix} at Wikipedia.
In this case the design matrix has two columns: the first has the x-values, the
second contains just ones.

@examples[#:eval quick-eval
                (define xs (column 0 1 2 3))
                (define ys (column 1 3 5 7))
                (define X  (augment xs (flmatrix-ones (nrows xs) 1)))
                X
                (define B  (lstsq X ys))
                B]

@section{Installation}

The package @racket[flmatrix] is installed either in the terminal:

    @racket[raco pkg install flmatrix]

or using the Package Manager in DrRacket.

The package relies on the shared libraries CBLAS and LAPACK.
Depending on your OS, you might need to install these yourself.

On macOS both CBLAS and LAPACK is part of the Accelerate Framework
which is distributed by Apple. This means no extra installation is
needed.

On Linux you need copies of CBLAS and LAPACK. Since BLAS and LAPACK
exists in multiple versions, so a little care is needed. First
on most systems @racket[libblas] is used for the Fortran version,
and @racket[libcblas], so get the latter. However on Debian it turns
out @racket[libblas] is exporting the names used by CBLAS, so
(either?) ought to be fine.

On Windows: A tester is needed. Install CBLAS and LAPACK and let
me know if it works. Otherwise make an Issue at Github and we
will add the proper paths.



