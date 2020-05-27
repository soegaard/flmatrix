#lang scribble/manual
@(require scribble/example)
@(define reference.scrbl '(lib "scribblings/reference/reference.scrbl"))
@(define math.scrbl      '(lib "math/scribblings/math.scrbl"))

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
  @item{Level 1:  High level   - do what I mean}
  @item{Level 2:  Median level - do what I mean this way}
  @item{Level 3:  Low level    - do it using this underlying C-function}]

To use the library one can get by with level 1 operations, but if you understand
the underlying representation, you can improve your algorithms using 
level 2 operations. For those that really want to squeeze out the last bit of
performance we have made level 3 operations available as well.


@section{Quick Tutorial}

@(define quick-eval (let ([e (make-base-eval)]) (e '(require flmatrix/flmatrix)) e))

This section shows how to do simple matrix computations.

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
   return a column vector with @racket[num] numbers evenly spaced from @racket[start] to @racket[stop]

@bold[@racket[(linspace start stop num #f)]] @linebreak[]
   like @racket[(linspace start stop num)] but omit the last number

@examples[#:eval quick-eval
          (linspace 2 4 6)
          (linspace 2 4 6 #f)]










