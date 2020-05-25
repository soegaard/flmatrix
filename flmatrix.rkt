#lang racket
(provide (all-defined-out))

;;; TODO   
;;;        * Improve matrix-expt! (avoid allocation)

;;; NOTES
;;;        * Contracts will be added before release
;;;
;;;        * See tests at bottom for examples.

;;; FEEDBACK
;;;        * Where is CBLAS and LAPACK on your platform
;;;          (Windows and Linux)
;;;        * What are the libraries named?
;;;        * Do all tests evaluate to #t on your platform?

;;;        * Mail: jensaxel@soegaard.net

;;;
;;; PLATFORMS TESTED       
;;;        * OS X Mountain Lion (Working)

;;;
;;; IDEAS
;;;      Potential Improvements
;;;        * DONE Unsafe operations
;;;        * DONE add lda to the flmatrix structure
;;;        * DONE support shared submatrix without allocation
;;;        * DONE Improve equal?
;;;        * Use dgeequ before dgetrf (in matrix-lu!)
;;;        * Use an extra call with lwork=-1 in matrix-inverse!
;;;        * support different storage schemes 
;;;          http://www.netlib.org/lapack/lug/node121.html

;;; Useful routines to consider:

;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlazro.html
;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlaset.html
;;;   Constructs diagonal matrices. Use for flmatrix-identity
;;; * http://www.math.utah.edu/software/lapack/lapack-d/dlaswp.html
;;;   Row interchanges
;;; * http://www.math.utah.edu/software/lapack/lapack-d/drscl.html
;;;   Scale by 1/a with correct rounding


(require ffi/vector
         ffi/unsafe
         ffi/unsafe/define
         racket/flonum
         (for-syntax 
          racket/format
          racket/string
          ffi/unsafe))

;;;
;;; LIBRARIES
;;;

; CBLAS and LAPACK are used.
; The first two are C-based whereas LAPACK is Fortran based.

; Note: Use trailing _ in names exported by LAPACK (in order to work both on macOS and Linux).

;; Find placement of libraries.

(define-values (cblas-lib lapack-lib)
  (case (system-type)
    ; MACOS
    [(macosx)
     (define veclib-lib 
       ; OS X: Contains CBLAS both CATLAS. CATLAS is not used here.
       ; https://developer.apple.com/library/mac/#documentation/Accelerate/
       ;         Reference/BLAS_Ref/Reference/reference.html
       (ffi-lib "/System/Library/Frameworks/vecLib.framework/Versions/Current/vecLib"))     
     (define cblas-lib veclib-lib)
     (define lapack-lib
       (ffi-lib 
        (string-append
         "/System/Library/Frameworks/Accelerate.framework/"
         "Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK")))
     (values cblas-lib lapack-lib)]
    ; UNIX
    [(unix)
     (define cblas-lib  (ffi-lib "libblas"   '("3" #f)))  ; works on debian
     (define lapack-lib (ffi-lib "liblapack" '("3" #f)))  ; works on debian
     (values cblas-lib lapack-lib)]
    [(windows)
     ; tester needed
     (error 'tester-needed)]))

;;; Load libraries

(define-ffi-definer define-cblas  cblas-lib)
(define-ffi-definer define-lapack lapack-lib)

;;;
;;; REFERENCES
;;; 

; LAPACK Naming scheme:
; http://www.netlib.org/lapack/lug/node24.html

;;;
;;; CONFIGURATION
;;;

(define epsilon 1e-13) 
; If two flmatrices have the same size and
; the differences between two entries are
; smaller than epsilon, they are considered
; equal? . Furthermore if all entries are
; smaller than epsilon flmatrix-zero? 
; returns true.

(define current-max-flmatrix-print-size (make-parameter 100))
; For matrices with smaller size, all
; entries are printed. For larger matrices
; only the dimension is printed.


;;;
;;; REPRESENTATION
;;;

; BLAS/LAPACK represents matrices as one-dimensional arrays
; of numbers (S=single, D=double, X=complex or Z=double complex).
; This library uses arrays of doubles.

(define _flmatrix (_cpointer 'flmatrix))

; The array is wrapped in a struct, which besides
; a pointer to the array, holds the number of
; rows and columns. Future extension could be to
; allow different types of numbers, or perhaps
; choose specialized operations for triangular matrices.

(define (flmatrix-print A port mode)
  (define print (if mode write display))
  (print 
   (if (< (flmatrix-size A) (current-max-flmatrix-print-size))
       ; constructor style printing:
       (list 'flmatrix: ; (flmatrix-m A) (flmatrix-n A) 
             (flmatrix->lists A))
       ; omit actual elements
       (list 'flmatrix (flmatrix-m A) (flmatrix-n A) 
             "..."))
   port))



(define (flmatrix= A B [eps #f])
  (define-param (m n a lda) A)
  (define-param (r c b ldb) B)
  (and (= m r) (= n c)
       (for*/and ([j (in-range n)]
                  [i (in-range m)])
         (define aij (unsafe-ref a lda i j))
         (define bij (unsafe-ref b ldb i j))
         (if eps
             (fl<= (flabs (fl- aij bij)) eps)
             (fl= aij bij)))))

; m = rows, n = cols, a = mxn array of doubles
; lda = leading dimension of a (see below)
(struct flmatrix (m n a lda)
  #:methods gen:custom-write 
  [(define write-proc flmatrix-print)]
  #:methods gen:equal+hash
  [(define equal-proc 
     (λ (A B rec)
       (and (= (flmatrix-m A) (flmatrix-m B))
            (= (flmatrix-n A) (flmatrix-n B))
            (or (equal? (flmatrix-a A) (flmatrix-a B))
                (flmatrix= A B epsilon)))))
   (define hash-proc  
     ; TODO: Avoid allocation in hash-proc.
     (λ (A rec) 
       (define-param (m n) A)
       (rec (cons m (cons n (flmatrix->vector A))))))
   (define hash2-proc 
     (λ (A rec) 
       (define-param (m n) A)
       (rec (cons n (cons m (flmatrix->vector A))))))]) 

; convenient destructuring
(define-syntax (define-param stx)
  (syntax-case stx ()
    [(_ (m n) A)
     #'(begin
         (define A1 A) 
         (define m (flmatrix-m A1))
         (define n (flmatrix-n A1)))]
    [(_ (m n a) A)
    #'(begin
        (define A1 A) 
        (define m (flmatrix-m A1))
        (define n (flmatrix-n A1))
        (define a (flmatrix-a A1)))]
    [(_ (m n a lda) A)
    #'(begin
        (define A1 A) 
        (define m   (flmatrix-m A1))
        (define n   (flmatrix-n A1))
        (define a   (flmatrix-a A1))
        (define lda (flmatrix-lda A1)))]
    [_
     (syntax/loc stx (error "Wrong number of arguments"))]))

;;;
;;; MEMORY LAYOUT
;;;

; The entries are layed out in column major order.
; This means that the entries in a column are
; contigious. LAPACK needs this order.

;  a[0]   a[0   +lda]  a[0   + 2*lda] ... a[0+(n-1)*lda]
;  a[1]   a[1   +lda]
;  a[2]
;  ...    ...
; a[m-1]  a[m-1 +lda]  a[m01 + 2*lda] ... a[m-1+(n-1)*lda] 

; For most matrices lda=m.

; For a submatrix it is possible that lda is larger than m.
; See http://stackoverflow.com/q/5009513/23567 
; Example:
;   If   ma=10, na=12, a=<some adress>, lda=10,
;   then mb=7,  nb=2,  b=a+3+4*lda, ldb=10 (=lda)
; represent a 7x2 submatrix whose upper, lefter
; corner in A is (3,4) (indices are 0-based).

; The array index of the (i,j)th entry is:
(define-syntax-rule (index lda i j)
    (+ i (* j lda)))

(define (ptr-elm a lda i j)
  ; address of (i,j)th element
  (ptr-add a (index lda i j) _double))

(define (shared-submatrix! A i j r s)
  ; return rxs matrix with upper left corner (i,j)
  ; entries are shared with A
  ; TODO: consider garbage collection
  (define-param (m n a lda) A)
  (flmatrix r s (ptr-elm a lda i j) lda))

(define (flsubmatrix A m n i j)
  ; TODO: argument order not consistent with shared-submatrix!
  ; return a the mxn submatrix of with upper 
  ; left corner in (i,j)
  (copy-flmatrix (shared-submatrix! A i j m n)))


(define (ptr-row a i)
  ; pointer to beginning of row a
  (ptr-add a i _double))

(define (ptr-col a lda j)
  ; address of column j
  (ptr-add a (* j lda) _double))


;;;
;;; CHECKS
;;;

(define (check-flmatrix who A)
  (unless (flmatrix? A)
    (raise-type-error who "expected flmatrix" A)))

(define (check-same-dimensions A B who)
  (unless (flmatrix-same-dimensions? A B)
    (raise-argument-error who "expected two matrices of the same size" A B)))

(define (check-product-dimensions who A B [C #f])
  (unless (if (not C)
              (= (flmatrix-n A) (flmatrix-m B))
              (and (= (flmatrix-n A) (flmatrix-m B))
                   (= (flmatrix-m A) (flmatrix-m C))
                   (= (flmatrix-n B) (flmatrix-n C))))
    (raise-argument-error 
     who 
     (if C
         "expected three matrices with compatible dimensions"
         "expected two matrices with compatible dimensions")
     A B C)))

(define (check-matrix-vector-product-dimensions who A X Y transpose-A)
  (define-param (ma na) A)
  (define-param (mx nx) X)
  (define-param (my ny) Y)
  (unless (and (= ny nx 1)
               (if transpose-A
                 (and (= ma mx) (= na my))
                 (and (= na mx) (= ma my))))
    (raise-argument-error 
     who
     "expected same number of rows"
     (list (list ma na)
           (list mx nx)
           (list my ny)
           (list 'transpose-A transpose-A)))))

(define (check-legal-column who j A)
  (unless (< j (flmatrix-n A))
    (raise-argument-error 
     who "column index too large" j))
  (unless (<= 0 j)
    (raise-argument-error 
     who "column index must be non-negative")))

(define (check-legal-row who i A)
  (unless (< i (flmatrix-m A))
    (raise-argument-error 
     who "row index too large" i))
  (unless (<= 0 i)
    (raise-argument-error 
     who "row index must be non-negative")))

(define (check-square who A)
  (define-param (m n) A)
  (unless (= m n)
    (raise-argument-error 
     who "square matrix expected" A)))

(define (check-vector who v)
  (unless (vector? v) (raise-argument-error who "vector expected" v)))

(define (check-integer who x)
  (unless (integer? x) (raise-argument-error who "integer expected" x)))



;;;
;;; SIZE and DIMENSION
;;;

(define (flmatrix-size A)
  (check-flmatrix 'flmatrix-size A )
  (define-param (m n) A)
  (* m n))

(define (flmatrix-dimensions A)
  (check-flmatrix 'flmatrix-dimensions A)
  (define-param (m n) A)
  (values m n))

(define (flmatrix-same-dimensions? A B)
  (define-param (ma na) A)
  (define-param (mb nb) B)
  (and (= ma mb) (= na nb)))

(define (flmatrix-row-vector? A)
  (= 1 (flmatrix-m A)))

(define (flmatrix-column-vector? A)
  (= 1 (flmatrix-n A)))

;;;
;;; ALLOCATIONS and CONSTRUCTORS
;;; 

(define (alloc-flmatrix m n)
  (if (or (= m 0) (= n 0))
      #f ; ~ NULL
      (cast (malloc (* m n) _double 'atomic-interior)
            _pointer _flmatrix)))

(define (alloc-same-size-matrix A)
  (define-param (m n) A)
  (alloc-flmatrix m n))


(define-syntax (define-cblas* stx)
  (syntax-case stx ()
    [(def xname _x (c ...) body ...)
     (let ()       
       (define ((xname->name ctx xname) c)
         (datum->syntax 
          ctx
          (string->symbol
           (string-replace (~a xname) "x" (~a c) #:all? #f))))
       (define (c->_c c)
         (unless (symbol? c)
           (error (format "expected symbol, got: ~a" c)))
         (case c
           [(c) #'_double] ; TODO missing from ffi?
           [(z) #'_double] ; TODO
           [(d) #'_double]
           [(s) #'_float]
           [else (error "expected one of c, z, d, s")]))
       (with-syntax ([(name ...) 
                      (map (xname->name stx (syntax->datum #'xname))
                           (syntax->datum #'(c ...)))]
                     [(_c ...) (map c->_c (syntax->datum #'(c ...)))])
         #'(begin
             (define-cblas name 
               (let ([_x _c]) body ...))
             ...)))]))

(define-cblas* cblas_xcopy _x (s d c z)
  ; copy n elements from vector X to vector Y
  (_fun (n : _int)  
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _void))

#;(define-cblas cblas_dcopy
  ; copy n elements from vector X to vector Y
  (_fun (n : _int)  
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _void))

(define (unsafe-vector-copy! s a lda b)
  ; copy s elements from A into B
  ; element 0, lda, 2*lda, ... is copied
  (cblas_dcopy s a lda b 1))

(define (unsafe-matrix-copy! m n a lda b ldb)
  ; copy the mxn matrix A into B
  ; copy has upper left corner in (i,j)
  ; Note: use (ptr-elm b ldb i j) to
  ;       copy into a submatrix of b.
  (for ([j (in-range n)])
    (unsafe-vector-copy! 
     m (ptr-elm a lda 0 j) 1 
     (ptr-add b (* j ldb) _double))))

(define (copy-flmatrix A)
  (define-param (m n a lda) A)
  (define size (* m n))
  (define b (cast (malloc size _double 'atomic)
                  _pointer _flmatrix))
  (define ldb m)
  (cond 
    [(= lda m) ; elements in a are contigious
     (unsafe-vector-copy! size a 1 b)]
    [else ; copy each column separately
     (unsafe-matrix-copy! m n a lda b ldb)])
  (flmatrix m n b ldb))

(define (make-flmatrix m n [x 0.0])
  (define a (alloc-flmatrix m n))
  (define x* (real->double-flonum x))
  (if (= x 0.0)
      (memset a 0 (* m n) _double)
      (for ([i (* m n)]) (ptr-set! a _double i x*)))
  (flmatrix m n a m))

(define (flmatrix-zeros m n)
  (make-flmatrix m n 0.0))

(define (flmatrix-ones m n)
  (make-flmatrix m n 1.0))


(define (list->flmatrix xss)
  (define m (length xss))
  (define n (apply max (map length xss)))
  (for*/flmatrix m n
                 ([xs (in-list xss)]
                  [x  (in-list xs)])
                 x))

(define (vectors->flmatrix xss)
  (define m (vector-length xss))
  (define n (vector-length (vector-ref xss 0)))
  (for*/flmatrix m n
                 ([xs (in-vector xss)]
                  [x  (in-vector xs)])
                 x))

(define (flmatrix-identity m)
  (define A (make-flmatrix m m 0.0))
  (for ([i (in-range m)])
    (flmatrix-set! A i i 1.0))
  A)

(define (flmatrix-column A j)
  ; copy column j
  (check-legal-column 'flmatrix-column j A)  
  (define-param (m n) A)
  (copy-flmatrix (shared-submatrix! A 0 j m 1)))

(define (flmatrix-row A i)
  ; copy row i
  (define-param (m n) A)
  (check-legal-row 'flmatrix-row i A)
  (copy-flmatrix (shared-submatrix! A i 0 1 n)))

;;;
;;; CONVERSIONS MATRIX <-> VECTOR
;;;

(define (flmatrix->vector A)
  ; the result vector uses row-major order
  (define-param (m n a lda) A)
  (for*/vector #:length (* m n)
    ([i (in-range 0 m)]
     [j (in-range 0 n)])
    (unsafe-ref a lda i j)))

(define (flmatrix->vectors A)
  ; the result is a vector of rows
  (define-param (m n a lda) A)
  (for/vector #:length m
    ([i (in-range 0 m)])
    (for/vector #:length n
        ([j (in-range 0 n)])
      (ptr-ref a _double (+ i (* j lda)))
      #;(ptr-ref (ptr-elm a lda i j) _double))))

(define (vector->flmatrix m n v)
  (unless (= (* m n) (vector-length v))
    (raise-argument-error
     'vector->flmatrix
     "expected m*n to be the same as the length of the vector"))
  (define a (alloc-flmatrix m n))
  (define k 0)
  (for* ([j (in-range n)]
         [i (in-range m)])
    (ptr-set! a _double* k ; (index m i j) 
              (vector-ref v (+ (* i n) j)))
    (set! k (+ k 1)))  
  (flmatrix m n a m))

; (: matrix/dim : Integer Integer Number * -> (Matrix Number))
; construct a mxn flmatrix with elements from the values xs
; the length of xs must be m*n
(define (flmatrix/dim m n . xs)
  (vector->flmatrix m n (list->vector xs)))


;;;
;;; COMPREHENSIONS
;;;

; (for/flmatrix m n (clause ...) . defs+exprs)
; (for/matrix (i in m) (j in n) (clauses ...) . body)
;    Return an  m x n  flmatrix with elements from the last expr.
;    The first n values produced becomes the first row.
;    The next n values becomes the second row and so on.
;    The bindings in clauses run in parallel.
(define-syntax (for/flmatrix stx)
  (syntax-case stx (in)
    [(_for/matrix (i in m-expr)  (j in n-expr) #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define a (alloc-flmatrix m n))
         (define idx 0)
         (for* ([j (in-range n)]
                [i (in-range m)]                
                clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx 1)))
         (flmatrix m n a m)))]
    ; elements in column 0 are generated first, then column 1, ...
    [(_ m-expr n-expr #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n  n-expr])
         (define a (alloc-flmatrix m n))
         (define size (* m n))
         (for ([idx (in-range size)] clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x))
         (flmatrix m n a m)))]
    [(_for/matrix (i in m-expr) (j in n-expr) (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define size (* m n))
         (define a (alloc-flmatrix m n))
         (define idx 0)
         (for* ([i (in-range m)]
                [j (in-range n)]                                
                clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flmatrix m n a m)))]
    ; elements in row 0 are generated first, then row 1, ...
    [(_ m-expr n-expr (clause ...) . defs+exprs)
     (syntax/loc stx
       (let* ([m m-expr] [n n-expr])
         (define a (alloc-flmatrix m n))
         (define idx 0)
         (define size (* m n))
         (for ([k (in-range size)] clause ...)
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flmatrix m n a m)))]))

; (for*/flmatrix m n (clause ...) . defs+exprs)
;    Return an  m x n  flmatrix with elements from the last expr.
;    The first n values produced becomes the first row.
;    The next n values becomes the second row and so on.
;    The bindings in clauses run nested.
; (for*/flmatrix m n #:column (clause ...) . defs+exprs)
;    Return an  m x n  flmatrix with elements from the last expr.
;    The first m values produced becomes the first column.
;    The next m values becomes the second column and so on.
;    The bindings in clauses run nested.

(define-syntax (for*/flmatrix stx)
  (syntax-case stx ()
    [(_ m-expr n-expr #:column (clause ...) . defs+exprs)
     (syntax/loc stx
       (let* ([m  m-expr] [n  n-expr])
         (define a (alloc-flmatrix m n))
         (define idx 0)
         (define size (* m n))
         (for* (clause ... #:break (= idx size))
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx 1)))
         (flmatrix m n a m)))]
    [(_ m-expr n-expr (clause ...) . defs+exprs)
     (syntax/loc stx
       (let ([m m-expr] [n n-expr])
         (define a (alloc-flmatrix m n))
         (define idx 0)
         (define size (* m n))
         (for* (clause ... #:final (= idx (- size 1)))
           (define x (let () . defs+exprs))
           (ptr-set! a _double* idx x)
           (set! idx (+ idx m))
           (when (>= idx size)
             (set! idx (+ idx 1 (- size)))))
         (flmatrix m n a m)))]))

(define-syntax (for/flmatrix-sum stx)
  (syntax-case stx ()
    [(_ (for:-clause ...) . defs+exprs)
     (syntax/loc stx
       (let ()
         (define sum #f)
         (for (for:-clause ...)
           (define a (let () . defs+exprs))
           (set! sum (if sum (flmatrix+ sum a) a)))
         sum))]))

;;;
;;; BINARY MATRIX OPERATIONS
;;; 

;;; MATRIX SUM AND DIFFERENCE

(define-cblas* cblas_xaxpy _x (s d #;c #;z)
  ; Y := αX+Y  ; X and Y are vectors
  ; If incX=3 then every 3rd element of X is used.
  (_fun (n : _int) (alpha : _x)
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _void))

#;(define-cblas cblas_daxpy
  ; Y := αX+Y  ; X and Y are vectors
  ; If incX=3 then every 3rd element of X is used.
  (_fun (n : _int) (alpha : _double) 
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _void))

(define (unsafe-vector-clear n a [lda 1])
  (cblas_daxpy n -1.0 a lda a lda))

; TODO: Allow adding row to different matrix!

(define (flmatrix-add-scaled-row! A i1 s i2)
  ; scale row i2 and add to row i1
  (check-legal-row 'matrix-add-scaled-row! i1 A)
  (check-legal-row 'matrix-add-scaled-row! i2 A)
  (define-param (m n a lda) A)
  (define rowi1 (ptr-row a i1))
  (define rowi2 (ptr-row a i2))
  (define s* (real->double-flonum s))
  (cblas_daxpy n s* rowi2 lda rowi1 lda)
  A)

(define (flmatrix-add-scaled-row A i1 s i2)
  (define B (copy-flmatrix A))
  (flmatrix-add-scaled-row! B i1 s i2)
  B)

(define (flmatrix-add-scaled-column! A j1 s j2)    
  (check-legal-row 'flmatrix-add-scaled-column! j1 A)
  (check-legal-row 'flmatrix-add-scaled-column! j2 A)
  (define-param (m n a lda) A)
  (define colj1 (ptr-col a lda j1))
  (define colj2 (ptr-col a lda j2))
  (define s* (real->double-flonum s))
  (cblas_daxpy m s* colj1 1 colj2 1)
  A)

(define (flmatrix-add-scaled-column A i1 s i2)
  (define B (copy-flmatrix A))
  (flmatrix-add-scaled-column! B i1 s i2)
  B)

(define (constant*flmatrix+flmatrix! alpha A B)
  ; B := αA+B  
  (define-param (m n a lda) A)
  (define-param (r s b ldb) B)
  (for ([j (in-range n)])
    (cblas_daxpy m alpha 
                 (ptr-col a lda j) 1 
                 (ptr-col b ldb j) 1))
  B)

(define (constant*flmatrix+flmatrix alpha A B)
  ; αA+B
  (define αA+B (copy-flmatrix B))
  (constant*flmatrix+flmatrix! alpha A αA+B)
  αA+B)

(define (flmatrix+! A B)
  ; B := A + B
  (check-same-dimensions A B 'flmatrix+!)
  (constant*flmatrix+flmatrix! 1.0 A B))

(define (flmatrix+ A B)
  ; A + B
  (check-same-dimensions A B 'flmatrix+)
  (constant*flmatrix+flmatrix 1.0 A B))

(define (flmatrix-! A B)
  ; A := A - B
  (check-same-dimensions A B 'flmatrix-!)
  (constant*flmatrix+flmatrix! -1.0 B A))

(define (flmatrix- A [B #f])
  (cond
    [B
     (check-same-dimensions A B 'flmatrix-)
     (constant*flmatrix+flmatrix -1.0 B A)]
    [else
     (flmatrix-scale -1.0 A)]))

;;; Matrix x Matrix Multiplication

(define _CBLAS_ORDER _int)
(define CblasRowMajor 101)
(define CblasColMajor 102)

(define _CBLAS_TRANSPOSE _int)
(define CblasNoTrans   111)
(define CblasTrans     112)
(define CblasConjTrans 113)

(define-cblas* cblas_xgemm _x (s d z c)
  ; C := α(A*B)+βC 
  ; 1. Multiplies A and B.
  ; 2. Scales result with alpha
  ; 3. Scales C with beta.
  ; 4. Stores sum in in C.
  (_fun (order : _CBLAS_ORDER) 
        (transa : _CBLAS_TRANSPOSE) ; transpose A?
        (transb : _CBLAS_TRANSPOSE) ; transpose B?
        (m : _int) ; rows in A and C
        (n : _int) ; cols in B and C
        (k : _int) ; cols in A = rows in B
        (alpha : _x) ; scaling factor for A and B
        (A : _flmatrix) 
        (lda : _int) ; size of first dim of A
        (B : _flmatrix) 
        (ldb : _int) ; size of first dim of B
        (beta : _double) ; scaling for C
        (C : _flmatrix) 
        (ldc : _int) ; size of first dim of C
        -> _void))

(define (constant*matrix*matrix+constant*matrix! alpha A B beta C transA transB)
  ; C := α(A*B)+βC, maybe transpose A and/or B first 
  (check-product-dimensions 'constant*matrix*matrix+constant*matrix! A B C)
  (define-param (m n a lda) A)
  (define-param (r s b ldb) B)
  (define-param (x y c ldc) C)
  (define alpha* (real->double-flonum alpha))
  (define beta*  (real->double-flonum beta))
  (cblas_dgemm CblasColMajor 
               (if transA CblasTrans CblasNoTrans)
               (if transB CblasTrans CblasNoTrans)
               m s n alpha* 
               a lda  b ldb  beta*  c ldc)
  C)

(define (flmatrix*! A B C 
                    [alpha 1.0] [beta 1.0] 
                    [transpose-A #f] [transpose-B #f])
  ; C := α(A*B)+βC, maybe transpose A and/or B first 
  (constant*matrix*matrix+constant*matrix! 
   alpha A B beta C transpose-A transpose-B))

(define (flmatrix* A B [C #f]
                   [alpha 1.0] [beta 1.0] 
                   [transpose-A #f] [transpose-B #f])
  ; C := α(A*B)+βC, maybe transpose A and/or B first 
  (define C1 (or C (make-flmatrix (flmatrix-m A) (flmatrix-n B))))
  (flmatrix*! A B C1 alpha beta transpose-A transpose-B))

;;; Matrix Power

(define (flmatrix-expt a n)
  (check-flmatrix 'flmatrix-expt a)
  (check-square 'matrix-expt a)
  (cond
    [(= n 0)  (flmatrix-identity (flmatrix-m a))]
    [(= n 1)  (copy-flmatrix a)]
    [(= n 2)  (flmatrix* a a)]
    [(even? n) (let ([a^n/2 (flmatrix-expt a (quotient n 2))])
                 (flmatrix* a^n/2 a^n/2))]
    [else     (flmatrix* a (flmatrix-expt a (sub1 n)))]))

;;; Matrix x Vector Multiplication

; NOTE: Functions accepting column vectors automatically
;       convert (standard) vectors into mx1 matrices.

(define-cblas* cblas_xgemv _x (s d c z) ; Double GEneral Matrix Vector multiplication
  ; Y := α(AX) +(βY) 
  (_fun (order : _CBLAS_ORDER) 
        (transa : _CBLAS_TRANSPOSE) ; transpose A?
        (m : _int) ; rows in A 
        (n : _int) ; cols in A 
        (alpha : _x) ; scaling factor for A 
        (A : _flmatrix) 
        (lda : _int) 
        (X : _flmatrix) ; vector
        (ldx : _int) 
        (beta : _x) ; scaling for Y
        (Y : _flmatrix) ; vector
        (ldy : _int) 
        -> _void))

(define (constant*matrix*vector+constant*vector! alpha A X beta Y transA)
  ; unsafe:  Y := α(AX) +(βY), maybe transpose A first 
  (define-param (m n a lda) A)
  (cblas_dgemv CblasColMajor 
               (if transA CblasTrans CblasNoTrans)
               m n
               (real->double-flonum alpha)
               a lda
               (flmatrix-a X) 1 
               (real->double-flonum beta)
               (flmatrix-a Y) 1)
  Y)


(define (flmatrix*vector! A X Y [alpha 1.0] [beta 1.0] 
                          [transpose-A #f])
  (define X1 (result-flcolumn X))
  (define Y1 (result-flcolumn Y))  
  (check-matrix-vector-product-dimensions 
   'constant*matrix*vector+constant*vector! A X1 Y1 transpose-A)
  ; Y := α(AX) +(βY), maybe transpose A first 
  (constant*matrix*vector+constant*vector! 
   alpha A X1 beta Y1 transpose-A))

(define (flmatrix*vector A X [Y #f] [alpha 1.0] [beta 1.0] 
                         [transpose-A #f] )
  ; Y := α(AX) +(βY), maybe transpose A first
  (define Y1 (or Y (make-flmatrix (if transpose-A (flmatrix-n A) (flmatrix-m A)) 1 0.0)))
  (flmatrix*vector! A X Y1 alpha 1.0 transpose-A))

;;;
;;; ELEMENT WISE OPERATIONS
;;;

;;; Ref

(define (unsafe-ref a lda i j)
  (ptr-ref (ptr-elm a lda i j) _double))

(define (flmatrix-ref A i j)
  (define-param (m n a lda) A)
  (unless (< -1 i m)
    (raise-arguments-error 
     'matrix-ref (format "expected row index between 0 and ~a, got ~a" m i)))
  (unless (< -1 j n)
    (error 'matrix-ref 
           (format "expected column index between 0 and ~a, got ~a" n j)))
  (unsafe-ref a lda i j))

;;; Set!

(define (unsafe-set! a lda i j x)
  (ptr-set! (ptr-elm a lda i j) _double x))

(define (flmatrix-set! A i j x)
  (check-legal-row    'flmatrix-set! i A)
  (check-legal-column 'flmatrix-set! j A)
  (define-param (m n a lda) A)
  (define x* (real->double-flonum x))
  (unsafe-set! a lda i j x*)
  A)

;;; Scaling

(define-cblas* cblas_xscal _x (s d c z)
  ; X := αX  vector
  (_fun (n : _int) (alpha : _x) 
        (X : _flmatrix) (incX : _int)        
        -> _void))

(define (constant*matrix! s A)
  ; A := s*A
  (define-param (m n a lda) A)
  (define s* (real->double-flonum s))
  (cond 
    [(= lda m)
     (cblas_dscal (* m n) s* a 1)]
    [else 
     (for ([j (in-range n)])
       (cblas_dscal m s* (ptr-col a lda j) 1))])
  A)

(define (flmatrix-scale! s A)
  ; A := s*A
  (constant*matrix! s A))

(define (flmatrix-scale s A)
  ; s*A
  (define sA (copy-flmatrix A))
  (flmatrix-scale! s sA))

(define (shared-column-flmatrix A j)
  (check-legal-column 'shared-column-flmatrix j A)
  (define-param (m n) A)
  (shared-submatrix! A 0 j m 1))

(define (shared-row-flmatrix A i)
  (check-legal-row 'shared-row-flmatrix i A)
  (shared-submatrix! A i 0 1 (flmatrix-n A)))

(define (flmatrix-scale-column! A j s)
  ; col_j := s * col_j
  (constant*matrix! s (shared-column-flmatrix A j))
  A)

(define (flmatrix-scale-column A j s)
  (define B (copy-flmatrix A))
  (flmatrix-scale-column! B j s)
  B)

(define (flmatrix-scale-row! A i s)
  ; row_i := s * rwo_i
  (check-legal-row 'flmatrix-scale-row! i A)
  (define-values (m n) (flmatrix-dimensions A))
  (constant*matrix! s (shared-row-flmatrix A i))
  A)

(define (flmatrix-scale-row A i s)
  (define B (copy-flmatrix A))
  (flmatrix-scale-row! B i s)
  B)

;;; Swapping

(define-cblas* cblas_xswap _x (s d c z)
  ; Swaps elements in the vectors x and y
  (_fun (n : _int) ; length of vector
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _void))

(define (flmatrix-swap-rows! A i1 i2)
  (check-legal-row 'flmatrix-swap-rows! i1 A)
  (check-legal-row 'flmatrix-swap-rows! i2 A)
  (unless (= i1 i2)
    (define-param (m n a lda) A)
    (define rowi1 (ptr-row a i1))
    (define rowi2 (ptr-row a i2))
    (cblas_dswap n rowi1 lda rowi2 lda))
  A)

(define (flmatrix-swap-rows A i1 i2)
  (define B (copy-flmatrix A))
  (flmatrix-swap-rows! B i1 i2)
  B)

(define (flmatrix-swap-columns! A j1 j2)
  (check-legal-row 'flmatrix-swap-columns! j1 A)
  (check-legal-row 'flmatrix-swap-columns! j2 A)
  (unless (= j1 j2)
    (define-param (m n a lda) A)
    (define colj1 (ptr-col a lda j1))
    (define colj2 (ptr-col a lda j2))
    (cblas_dswap m colj1 1 colj2 1))
  A)

(define (flmatrix-swap-columns A j1 j2)
  (define B (copy-flmatrix A))
  (flmatrix-swap-columns! B j1 j2)
  B)

(define (flmatrix-flip-left-to-right! A)
  (check-flmatrix 'flmatrix-flip-left-to-right! A)
  (define-param (m n a lda) A)
  (unless (= n 1)
    (for ([j (in-range (quotient n 2))])
      (flmatrix-swap-columns! A j (- n j 1))))
  A)

(define (flmatrix-flip-left-to-right A)
  (check-flmatrix 'flmatrix-flip-left-to-right A)
  (define B (copy-flmatrix A))
  (flmatrix-flip-left-to-right! B)
  B)

(define (flmatrix-flip-up-to-down! A)
  (check-flmatrix 'flmatrix-flip-up-to-down! A)
  (define-param (m n a lda) A)
  (unless (= m 1)
    (for ([i (in-range (quotient m 2))])
      (flmatrix-swap-rows! A i (- m i 1))))
  A)

(define (flmatrix-flip-up-to-down A)
  (check-flmatrix 'flmatrix-flip-up-to-down A)
  (define B (copy-flmatrix A))
  (flmatrix-flip-up-to-down! B)
  B)

(define (flmatrix-rotate90 A)
  ; 1 2 3          3 6 9
  ; 4 5 6 becomes  2 5 8
  ; 7 8 9          1 4 7
  (check-flmatrix 'flmatrix-rotate90 A)
  (check-square   'flmatrix-rotate90 A)
  (define-param (m n a lda) A)  
  (for*/flmatrix n m
                 ([i (in-range (- n 1) -1 -1)]
                  [j (in-range (- m 1) -1 -1)])
                 (flmatrix-ref A (- n j 1) i)))




;;; Max Absolute Value

(define-cblas* cblas_ixamax _x (s d c z)
  ; Returns the index of the element with the largest 
  ; absolute value in a vector.
  (_fun (n : _int) (X : _flmatrix) (incX : _int)
        -> _int))


(define (flmatrix-max-abs-index A)
  (define-param (m n a lda) A)
  (cond
    [(= m lda) 
     (define idx (cblas_idamax (* m n) a lda))
     (values (remainder idx m) (quotient idx m))]
    [(= n 1)
     (define idx (cblas_idamax m a 1))
     (values (- idx 1) 0)]
    [else
     (define idx (make-vector n))
     (for ([j (in-range n)])
       (define i (cblas_idamax m (ptr-col a lda j) 1))
       (vector-set! idx j (cons (cons i j) (unsafe-ref a lda i j))))
     (define ij (car (vector-argmax cdr idx)))
     (values (car ij) (cdr ij))]))

(define (flmatrix-max-abs-value A)
  (define-values (i j) (flmatrix-max-abs-index A))
  (flmatrix-ref A i j))

(define (flmatrix-zero? A [eps epsilon])
  ; set eps=#f to use normal equal?
  (define val (flmatrix-max-abs-value A))
  (if eps (< (abs val) eps) (zero? val)))

(define (flmatrix-ones? A [eps #f])
  ; is A a square matrix with ones on the main diaginal and zero elsewhere?
  (define-param (m n a lda) A)  
  (and (= m n)
       (for*/and ([j (in-range n)]
                  [i (in-range m)])
         (define aij (unsafe-ref a lda i j))
         (if (= i j)
             (if eps
                 (fl<= (flabs (fl- aij 1.0)) eps)
                 (fl= aij 1.0))
             (if eps
                 (fl<= (flabs aij) eps)
                 (fl= aij 0.0))))))



;;;
;;; BLOCK LEVEL OPERATIONS
;;;

(define (flmatrix-augment C . Cs)
  ; 1. Check that all have same number of rows.
  (define-param (mc nc c ldc) C)
  (define rows (map flmatrix-m (cons C Cs)))
  (unless (andmap (λ (r) (= mc r)) rows)
    (raise-arguments-error 
     'flmatrix-augment
     "all arguments must have same number of rows"))
  ; 2. Find size for result matrix and allocate
  (define m mc)
  (define n (apply + (map flmatrix-n (cons C Cs))))
  (define a (alloc-flmatrix m n))
  (define lda m)
  ; 3. Fill in blocks
  (define j 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define aj (ptr-col a lda j))
    (unsafe-matrix-copy! mb nb b ldb aj lda)
    (set! j (+ j nb)))
  (flmatrix m n a lda))

(define (flmatrix-stack C . Cs)
  ; 1. Check that all have same number of columns
  (define-param (mc nc c ldc) C)
  (define cols (map flmatrix-n (cons C Cs)))
  (unless (andmap (λ (x) (= x nc)) cols)
    (raise-arguments-error 
     'flmatrix-stack
     "all arguments must have same number of columns"))
  ; 2. Find size for result matrix and allocate
  (define rows (map flmatrix-m (cons C Cs)))
  (define m (apply + rows))
  (define n nc)
  (define a (alloc-flmatrix m n))
  (define lda m)  
  ; 3. Fill in blocks
  (define i 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define ai (ptr-row a i))
    (unsafe-matrix-copy! mb nb b ldb ai lda)
    (set! i (+ i mb)))
  (flmatrix m n a lda))

(define (flmatrix-block-diagonal C . Cs)
  (define rows (map flmatrix-m (cons C Cs)))
  (define cols (map flmatrix-n (cons C Cs)))
  ; 2. Find size for result matrix and allocate
  (define m (apply + rows))
  (define n (apply + cols))
  (define a (alloc-flmatrix m n))
  (define lda m)
  (unsafe-vector-clear (* m n) a)
  ; 3. Fill in blocks
  (define i 0)
  (define j 0) 
  (for ([B (in-list (cons C Cs))])
    (define-param (mb nb b ldb) B)
    (define aij (ptr-elm a lda i j))
    (unsafe-matrix-copy! mb nb b ldb aij lda)
    (set! i (+ i mb))
    (set! j (+ j nb)))
  (flmatrix m n a lda))

(define (flmatrix-repeat A m [n m])
  ; Make a matrix with mxn blocks, each block is A.
  (define row (apply flmatrix-augment (for/list ([i m]) A)))
  (apply flmatrix-stack (for/list ([j n]) row)))
  
  

;;;
;;; NORMS
;;;

(define-cblas* cblas_xnrm2 _x (s d)
  ; L2-norm = (sqrt (sum (sqr X_i))), vector
  (_fun (n : _int) (X : _flmatrix) (incX : _int)
        -> _x))

(define (flmatrix-norm A)
  ; (sqrt (sum (sqr A_ij)))
  (define-param (m n a lda) A)
  (cond
    [(= lda m)
     (cblas_dnrm2 (* m n) a 1)]
    [(= n 1)
     (cblas_dnrm2 m a 1)]
    [else
     (sqrt
      (for/sum ([j (in-range n)])
        (expt (cblas_dnrm2 m (ptr-col a lda j) 1) 2)))]))

;;; 
;;; UNARY MATRIX OPERATIONS
;;;

(define (flmatrix-transpose A)
  ; TODO: Measure: Is it faster to use 
  ;       a loop with unsafe-vector-copy ?
  (define-param (m n a lda) A)
  (define AT (make-flmatrix n m))
  (define at (flmatrix-a AT))
  (for* ([j (in-range n)]
         [i (in-range m)])
    (unsafe-set! at n j i (unsafe-ref a lda i j)))
  AT)

;;;
;;; MATRIX DECOMPOSITIONS
;;;

;;; Pivots 

(struct pivots (ps)) ; ps is a u32vector
;   ps[i]=j  <=>  row i and row j-1 is swapped
; Note: Fortran counts from 1 !

(define (unsafe-pivot-ref ps i)
  ; Fortran indices are 1-based.
  (- (u32vector-ref ps i) 1))

(define (pivots-ref Ps i)
  (unsafe-pivot-ref (pivots-ps Ps) i))

(define (pivots-length Ps)
  (u32vector-length (pivots-ps Ps)))

(define (pivots->flmatrix Ps)
  ; return the permuation matrix
  (define ps (pivots-ps Ps))
  (define k (u32vector-length ps))  
  (define A (make-flmatrix k k 0.0))
  (define-param (m n a lda) A)
  ; introduce ones on diagonal
  (for ([i (in-range m)])
    (unsafe-set! a lda i i 1.0))
  ; perform row permutations
  (for ([i (in-range (- m 1) -1 -1)])
    (define i* (unsafe-pivot-ref ps i))
    (unless (= i i*)
      (flmatrix-swap-rows! A i i*)))
  A)

(define (pivots-sign Ps)
  ; return the sign of the corresponding permuation
  (define ps (pivots-ps Ps))
  (define n (u32vector-length ps))
  (for/product ([i (in-range n)])
    (define i* (unsafe-pivot-ref ps i))
    (if (= i i*) 1 -1)))

;;;
;;; PLU Factorization
;;;

; A = P L U
; where P is a permutation matrix,
;       L is lower triangular
;       U is upper triangular.
; Note: U is the result of Gauss elimation.

(define-lapack dgetrf_ 
  ; http://www.netlib.org/lapack/double/dgetrf.f
  ; DGETRF computes an LU factorization of a general M-by-N matrix A
  ; using partial pivoting with row interchanges.
  ; The factorization has the form
  ;     A = P * L * U
  ; where P is a permutation matrix, L is lower triangular with unit
  ; diagonal elements (lower trapezoidal if m > n), and U is upper
  ; triangular (upper trapezoidal if m < n).
  
  ; Algorithm: Gaussian elimination with partial pivoting
  (_fun (m : (_ptr i _int))
        (n : (_ptr i _int))
        (a : _flmatrix)
        (lda : (_ptr i _int))
        (ipiv : (_u32vector o (ptr-ref m _int)))
        (info : (_ptr o _int)) 
        -> _void
        -> (values (pivots ipiv) info)))

(define (flmatrix-lu! A)
  (define-param (m n a lda) A)
  (dgetrf_ m n a lda))

(define (flmatrix-plu A)
  (define B (copy-flmatrix A))
  (define-values (ps info) (flmatrix-lu! B))
  (define P (pivots->flmatrix ps))
  (define L (flmatrix-extract-lower B))
  (define U (flmatrix-extract-upper B))
  ; TODO: What to do with info?
  (values P L U))

(define (flmatrix-extract-upper A)
  ; extract the upper matrix, 
  ; including the diagonal
  ; discard below diagonal
  (define-param (m n) A)
  (define k (min m n))
  (define U (make-flmatrix k k))
  ; TODO: use unsafe-ref or unsafe-vector-copy
  (for* ([j (in-range 0 k)]
         [i (in-range (min (+ 1 j) k))])
    (flmatrix-set! U i j (flmatrix-ref A i j)))
  U)

(define (flmatrix-extract-lower A)
  ; extract the lower matrix, 
  ; and insert ones on diagonal
  (define L (copy-flmatrix A))
  (define-param (m n) A)
  ; TODO: use unsafe-ref or unsafe-vector-copy
  (for* ([j (in-range n)]
         [i (in-range 0 j)])
    (flmatrix-set! L i j 0))
  (for* ([j (in-range (min m n))])
    (flmatrix-set! L j j 1.0))
  L)

;;; SVD - Singular Value Decomposition

(define-lapack dgesvd_ 
  ; compute SVD 
  ; A = U * SIGMA * transpose(V)
  ; SIGMA is an mxm matrix, 
  ; Algorith: QR used
  (_fun (jobu : (_ptr i _byte))  ; char: a, s, o or n
        (jobvt : (_ptr i _byte)) ; char
        (m : (_ptr i _int)) ; rows in A
        (n : (_ptr i _int)) ; cols in A
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (s : _flmatrix) ; min(m,n) x 1
        (u : _flmatrix) ; mxm if jobu = a
        (ldu : (_ptr i _int)) 
        (vt : _flmatrix) ; nxn if jobvt = a
        (ldvt : (_ptr i _int)) ;         
        (work : _flmatrix) ; dim max(1,lwork)
        (lwork : (_ptr i _int)) ; 
        (info : (_ptr o _int))
        -> _void
        -> info))

(define-lapack dgesdd_ 
  ; compute SVD 
  ; A = U * SIGMA * transpose(V)
  ; SIGMA is an mxm matrix, 
  ; Algorithm: Divide and conquer with QR used for small
  ; This is the recommended algorithm, but uses
  ; more work space.
  (_fun (jobu : (_ptr i _byte)) ; char: a, s, o or n
        (jobvt : (_ptr i _byte))
        (m : (_ptr i _int)) ; rows in A
        (n : (_ptr i _int)) ; cols in A
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (s : _flmatrix) ; min(m,n) x 1
        (u : _flmatrix) ; mxm if jobu = a
        (ldu : (_ptr i _int)) 
        (vt : _flmatrix) ; nxn if jobvt = a
        (ldvt : (_ptr i _int)) ;         
        (work : _flmatrix) ; dim max(1,lwork)
        (lwork : (_ptr i _int)) ; 
        (info : (_ptr o _int))
        -> _void
        -> info))

(define (flmatrix-svd! A)
  ; TODO: Use lwork=-1 to get size of work
  (define-param (m n a lda) A)
  (define superb (- (min m n) 1))
  (define U  (make-flmatrix m m))
  (define S  (make-flmatrix (min m n) 1))
  (define VT (make-flmatrix n n))
  (define u  (flmatrix-a U))
  (define s  (flmatrix-a S))
  (define vt (flmatrix-a VT))
  (define lwork (* 10 (max m n))) ; conservative estimate
  (define W  (make-flmatrix lwork lwork))
  (define w  (flmatrix-a W))
  (define ca (char->integer #\A))
  (define info (dgesvd_ ca ca m n a lda s u m vt n w lwork))
  ; ? TODO: Best way to return error ?
  ; (when (> info 0) (displayln "Warning: no convergence"))
  (values U S VT))

(define (flmatrix-svd A)
  (flmatrix-svd! (copy-flmatrix A)))

(define (flmatrix-eigenvalues A)
  (define B (copy-flmatrix A))
  (define-values (S V D) (flmatrix-svd! B))
  (flmatrix->vector V))

(define (flmatrix-eigenvectors A)
  (define B (copy-flmatrix A))
  (define-values (S V D) (flmatrix-svd! B))
  (flmatrix->vector V))



;;; QR Factorization
; dgeqrfp returns positive entries on the diagonal
; for some reason this is missing on macOS, so now dgeqrf is used instead
#;(define-lapack dgeqrfp_ 
  ; Compute A = Q*R  
  ; Use dorgqr to generate matrix from output
  (_fun (m : (_ptr i _int)) ; rows in A
        (n : (_ptr i _int)) ; cols in A
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (tau : _flmatrix) ; min(m,n)x1        
        (work : _flmatrix) ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int)) ; >=max(1,n) best with >=n * blocksize
        (info : (_ptr o _int))  ; 
        -> _void
        -> info))

(define-lapack dgeqrf_
  ; Compute A = Q*R  
  ; Use dorgqr to generate matrix from output
  (_fun (m : (_ptr i _int)) ; rows in A
        (n : (_ptr i _int)) ; cols in A
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (tau : _flmatrix) ; min(m,n)x1        
        (work : _flmatrix) ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int)) ; >=max(1,n) best with >=n * blocksize
        (info : (_ptr o _int))  ; 
        -> _void
        -> info))


(define-lapack dorgqr_
  ; generate matrix from output of dgeqrf
  (_fun (m : (_ptr i _int)) ; rows in Q
        (n : (_ptr i _int)) ; cols in Q m>=n>=0
        (k : (_ptr i _int)) ; number of reflectors
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (tau : _flmatrix) ; min(m,n)x1        
        (work : _flmatrix) ; dim max(1,lwork) (x1)
        (lwork : (_ptr i _int)) ; >=max(1,n) best with >=n * blocksize
        (info : (_ptr o _int))  ; 
        -> _void
        -> info))


(define (flmatrix-qr B)
  (define A (copy-flmatrix B))
  (define-param (m n a lda) A)
  (define k (min m n))
  (define tau (make-flmatrix k k))
  (define atau (flmatrix-a tau))
  ; first call dgeqrf_ to get a working size
  (define work0 (make-flmatrix 1 1))
  (define info0 (dgeqrf_ m n a lda atau (flmatrix-a work0) -1))
  (define lwork (inexact->exact (flmatrix-ref work0 0 0))) ; 64 is a typical value
  ; now make the real call
  (define work  (make-flmatrix lwork 1))  
  (define awork (flmatrix-a work))  
  (define info (dgeqrf_ m n a lda atau awork lwork))
  (define R (flmatrix-extract-upper A))
  (define info1 (dorgqr_ m n k a lda atau awork lwork))  
  ; ? TODO: what to do with info  
  (values A R))

; old version used dgeqrfp
#;(define (flmatrix-qr B)
  (define A (copy-flmatrix B))
  (define-param (m n a lda) A)
  (define k (min m n))
  (define tau (make-flmatrix k k))
  (define atau (flmatrix-a tau))
  (define lwork (* 64 n)) ; 64 block size guess
  (define work (make-flmatrix lwork 1))
  (define awork (flmatrix-a work))
  ; TODO: Use lwork=-1 to get optimal lwork size
  (define info (dgeqrf_ m n a lda atau awork lwork))
  (define R (flmatrix-extract-upper A))
  (define info1 (dorgqr_ m n k a lda atau awork lwork))  
  ; ? TODO: what to do with info
  (values A R))

;;; 
;;; INVERSE
;;;

(define-lapack dgetri_
  ; http://www.netlib.org/lapack/double/dgetri.f
  ; DGETRI computes the inverse of a matrix using the LU factorization
  ; computed by DGETRF.
  ; This method inverts U and then computes inv(A) by solving the system
  ;   inv(A)*L = inv(U) for inv(A).
  (_fun (n : (_ptr i _int))
        (a : _flmatrix)
        (lda : (_ptr i _int))
        (ipiv : _u32vector)
        (work : (_or-null _flmatrix)) ; output
        (lwork : (_ptr i _int))
        (info : (_ptr o _int))
        -> _void
        -> (values info work)))

(define (flmatrix-inverse! A)
  ; TODO: this works, but call dgetri with lwork=-1
  ;       to get optimal size of workspace in first
  ;       entry of the work array.
  (define-param (m n a lda) A)
  (define work (copy-flmatrix A))
  (define-values (ipiv info) (flmatrix-lu! A))
  (dgetri_ m a lda (pivots-ps ipiv) (flmatrix-a work) (* m m))
  A)

(define (flmatrix-inverse A)
  (flmatrix-inverse! (copy-flmatrix A)))

;;;
;;; INVARIANTS
;;; 

(define (flmatrix-trace A)
  (check-square 'matrix-trace A)
  (for/sum ([i (in-range (flmatrix-m A))])
    (flmatrix-ref A i i)))

(define (flmatrix-determinant-from-plu LU pivots)
  ; compute determinant using output from PLU
  ; factorization
  (* (pivots-sign pivots)
     (for/product ([i (in-range (flmatrix-m LU))])
       (flmatrix-ref LU i i))))

(define (flmatrix-determinant A)
  (check-square 'matrix-determinant A)
  (define LU (copy-flmatrix A))
  (define-values (pivots info) (flmatrix-lu! LU))
  (flmatrix-determinant-from-plu LU pivots))

(define (count-columns-without-pivot pivots)
  ; TODO: Does this strategy work?
  (define ps (pivots-ps pivots))
  (define m (u32vector-length ps))
  (define with 
    (for/sum ([i (in-range m)])
      (define i* (- (u32vector-ref ps i) 1))
      (if (= i i*) 0 1)))
  (- m with))

(define (flmatrix-rank A)
  ; See answer: http://scicomp.stackexchange.com/questions/1861/understanding-how-numpy-does-svd
  ; rank = dimension of column space = dimension of row space  
  ;      = number of non-zero singular values
  (define-values (U Σ VT) (flmatrix-svd A))
  ; ? TODO: check info from -svd...
  (for/sum ([i (in-range (flmatrix-m Σ))]
            ; TODO: Which value for epsilon is correct?
            #:unless (< (abs (flmatrix-ref Σ i 0)) epsilon))
    1))

(define (flmatrix-nullity A)
  ; nullity = dimension of null space
  (define-param (m n) A)
  (- n (flmatrix-rank A)))

;;;
;;; VECTOR OPERATIONS
;;;

; Column vectors are represented as mx1 matrices.
; All operations working on column vectors accept
; standard vectors as input. Outputs are always
; in the form of a mx1 matrix.

(define (vector->flcolumn v)
  (define m (vector-length v))
  (vector->flmatrix m 1 v))

(define (vector->flrow v)
  (define n (vector-length v))
  (vector->flmatrix 1 n v))

(define (result-flcolumn c)
  ; convert output to mx1 matrix
  (if (vector? c)
      (vector->flcolumn c)
      c))

(define (flcolumn . xs)
  ; TODO: skip intermediary vector
  (vector->flcolumn 
   (list->vector xs)))

(define (flcolumn-size v)
  (if (vector? v)
      (vector-length v)
      (flmatrix-m v)))

;;; Dot Product

(define-cblas* cblas_xdot _x (s d)
  ; dot product, vectors
  (_fun (n : _int) 
        (X : _flmatrix) (incX : _int)
        (Y : _flmatrix) (incY : _int)
        -> _x))

(define (unsafe-vector-product n x y)
  (cblas_ddot n x 1 y 1))

(define (flcolumn-dot X Y)
  (set! X (result-flcolumn X))
  (set! Y (result-flcolumn Y))
  (define-param (m _ x ldx) X)
  (define-param (s __ y ldy) Y)
  (unless (= m s) 
    (error 
     'column-dot 
     "expected two mx1 matrices with same number of rows, got ~a and ~a"
     X Y))
  (unsafe-vector-product m x y))

(define fldot flcolumn-dot)

(define (flcolumn-norm v)
  (define-param (m _ a lda) (result-flcolumn v))
  (cblas_dnrm2 m a 1))

(define (flcolumn-unit m i)
  ; return i'th unit vector 
  (define U (make-flmatrix m 1 0.0))
  (flmatrix-set! U i 0 1.0)
  U)

(define (flscale-column s A)
  (define s* (real->double-flonum s))
  (cond
    [(vector? A)
     (define m (vector-length A))
     (vector->flcolumn 
      (for/vector #:length m 
        ([i (in-range m)])
        (* s* (vector-ref A i))))]
    [else
     (flmatrix-scale s A)]))

(define (flcolumn+ v w)
  (define m (flcolumn-size v))
  (define n (flcolumn-size w))
  (unless (= m n)
    (error 
     'flcolumn+ 
     "expected two column vectors of the same length, got ~a and ~a" v w))
  (cond 
    [(and (vector? v) (vector? w))
     (vector->flcolumn
      (for/vector #:length (+ m n) 
        ([i (in-range 0 m)]
         [x (in-vector v)]
         [y (in-vector w)])
        (+ x y)))]
    [else          
     (flmatrix+ (result-flcolumn v) (result-flcolumn w))]))

(define (flcolumn-projection v w)  
  ; Return the projection og vector v on vector w.
  (let ([w.w (fldot w w)])
    (if (zero? w.w)
        (error 'flcolumn-projection "projection on the zero vector not defined")
        (flscale-column (/ (fldot v w) w.w) w))))

(define (flcolumn-projection-on-unit v w)
  ; Return the projection of vector v on a unit vector w.  
  (flscale-column (flcolumn-dot v w) w))

(define (flcolumn-normalize w)
  ; Return unit vector with same direction as v.
  ; If v is the zero vector, the zero vector is returned.
  (define norm (flcolumn-norm w))
  (cond [(zero? norm) w]
        [else (flscale-column (/ norm) w)]))

(define (flzero-column-vector? v [eps #f])
  (define val (flmatrix-max-abs-value (result-flcolumn v)))
  (if eps (< (abs val) eps) (zero? val)))
  
; (flprojection-on-orthogonal-basis v bs)
;     Project the vector v on the orthogonal basis vectors in bs.
;     The basis bs must be either the column vectors of a matrix
;     or a sequence of column-vectors.
(define (flprojection-on-orthogonal-basis v bs)
  (if (empty? bs)
      (error 'flprojection-on-orthogonal-basis 
             "received empty list of basis vectors")
      (for/flmatrix-sum([b (in-list bs)])
                       (flcolumn-projection v (result-flcolumn b)))))

;     Project the vector v on the orthonormal basis vectors in bs.
;     The basis bs must be either the column vectors of a matrix
;     or a sequence of column-vectors.
(define (flprojection-on-orthonormal-basis v bs)
  (for/flmatrix-sum 
   ([b bs]) 
   (flmatrix-scale (flcolumn-dot v b) b)))
  

; (flgram-schmidt-orthogonal ws)
;     Given a list ws of flcolumn vectors, produce 
;     an orthogonal basis for the span of the
;     vectors in ws.
(define (flgram-schmidt-orthogonal ws1)
  (define ws (map result-flcolumn ws1))
  (cond 
    [(null? ws)       '()]
    [(null? (cdr ws)) (list (car ws))]
    [else 
     (define (loop vs ws)
       (cond [(null? ws) vs]
             [else
              (define w (car ws))
              (let ([w-proj (flprojection-on-orthogonal-basis w vs)])
                ; Note: We project onto vs (not on the original ws)
                ;       in order to get numerical stability.
                (let ([w-minus-proj (flmatrix- w w-proj)])
                  (if (flzero-column-vector? w-minus-proj)
                      (loop vs (cdr ws)) ; w in span{vs} => omit it
                      (loop (cons (flmatrix- w w-proj) vs) (cdr ws)))))]))
     (reverse (loop (list (car ws)) (cdr ws)))]))

; (flgram-schmidt-orthonormal ws)
;     Given a list ws of flcolumn vectors, produce 
;     an orthonormal basis for the span of the
;     vectors in ws.
(define (flgram-schmidt-orthonormal ws)
  (map flcolumn-normalize
       (flgram-schmidt-orthogonal ws)))

; (flprojection-on-subspace v ws)
;  Returns the projection of v on span{w_i}, w_i in ws.
(define (flprojection-on-subspace v ws)
  (flprojection-on-orthogonal-basis
   v (flgram-schmidt-orthogonal ws)))

;;;
;;; EQUATION SOLVING
;;;

(define-lapack dgesv_ ; Double, GEneral, Solve ...
  ; Compute solution to AX=B, where
  ; A is nxn and X and B are n x nrhs
  (_fun (n : (_ptr i _int))
        (nrhs : (_ptr i _int))
        (a : _flmatrix) ; io
        (lda : (_ptr i _int))
        (ipiv : (_u32vector o (ptr-ref n _int)))
        (b : _flmatrix) ; io
        (ldb : (_ptr i _int))
        (info : (_ptr o _int))
        -> _void
        -> info))

(define (flmatrix-solve A b)
  ; A matrix, b flcolumn
  (define-param (m n a lda) (copy-flmatrix A))
  (define bout (copy-flmatrix (result-flcolumn b)))
  (define info (dgesv_ n 1 a lda (flmatrix-a bout) m))
  ; ? TODO Handle info
  bout)

(define (flmatrix-solve-many! A B)  
  ; A matrix, b flcolumn  
  ; A and B are overwritten
  (define-param (m n a lda) A)
  (define-param (_ nrhs b ldb) B)
  (define info (dgesv_ n nrhs a lda b ldb))
  ; ? TODO: handle info
  (values B))

(define (flmatrix-solve-many A bs-or-B)  
  ; A matrix, b flcolumn 
  (define-param (m n) A)
  (define B (if (list? bs-or-B)    
                (apply flmatrix-augment 
                       (map result-flcolumn bs-or-B))
                (copy-flmatrix bs-or-B)))
  (flmatrix-solve-many! (copy-flmatrix A) B))

(define (flmatrix->columns A)
  (define-param (m n) A)  
  (for/list ([j (in-range n)])
    (flmatrix-column A j)))

;;;
;;; SEQUENCES
;;; 

(define (in-flrow/proc A r)
  (define-param (m n a lda) A)
  (make-do-sequence
   (λ ()
     (define (pos->elm j) (unsafe-ref a lda r j))
     (define next-pos add1)
     (define initial-pos 0)
     (define (continue? j) (< j n))
     (values pos->elm initial-pos continue? #f #f))))

; (in-flrow M i]
;     Returns a sequence of all elements of row i,
;     that is xi0, xi1, xi2, ...
(define-sequence-syntax in-flrow
  (λ () #'in-flrow/proc)
  (λ (stx)
    (syntax-case stx ()
      [[(x) (_ M-expr r-expr)]
       #'((x)
          (:do-in
           ([(M r m n a lda)
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 r-expr rd cd a lda))])
           (begin 
             (unless (flmatrix? M) 
               (raise-type-error 'in-flrow "expected flmatrix, got ~a" M))             
             (unless (and (integer? r) (and (<= 0 r ) (< r m))) 
               (raise-type-error 'in-flrow "expected row number" r)))
           ([j 0])
           (< j n)
           ([(x) (unsafe-ref a lda r j)])
           #true
           #true
           [(+ j 1)]))]
      [[(i x) (_ M-expr r-expr)]
       #'((i x)
          (:do-in
           ([(M r m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 r-expr rd cd a lda))])
           (begin 
             (unless (flmatrix? M) 
               (raise-type-error 'in-flrow "expected flmatrix, got ~a" M))             
             (unless (and (integer? r) (and (<= 0 r ) (< r m))) 
               (raise-type-error 'in-flrow "expected row number" r)))
           ([j 0])
           (< j n)
           ([(x) (unsafe-ref a lda r j)]
            [(i) j])
           #true
           #true
           [(+ j 1)]))]
      [[_ clause] (raise-syntax-error 
                   'in-flrow "expected (in-flrow <flmatrix> <row>)" #'clause #'clause)])))

; (in-flcol M j]
;     Returns a sequence of all elements of column j,
;     that is x0j, x1j, x2j, ...

(define (in-flcolumn/proc A s)
  (define-param (m n a lda) A)
  (make-do-sequence
   (λ ()
     (define (pos->elm i) (unsafe-ref a lda i s))
     (define next-pos add1)
     (define initial-pos 0)
     (define (continue? i) (< i m))
     (values pos->elm next-pos initial-pos #f #f))))

(define-sequence-syntax in-flcolumn
  (λ () #'in-flcolumn/proc)
  (λ (stx)
    (syntax-case stx ()
      ; M-expr evaluates to column
      [[(x) (_ M-expr)]
       #'((x)
          (:do-in
           ([(M n m a) 
             (let ([M1 (result-flcolumn M-expr)])
               (define-param (rd cd a) M1)
               (values M1 rd cd a))])
           (unless (flmatrix? M) 
             (raise-type-error 'in-column "expected matrix, got ~a" M))
           ([j 0])
           (< j n)
           ([(x) (ptr-ref a _double j)])
           #true
           #true
           [(+ j 1)]))]
      ; M-expr evaluates to matrix, s-expr to column index
      [[(x) (_ M-expr s-expr)]
       #'((x)
          (:do-in
           ([(M s m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 s-expr rd cd a lda))])
           (begin 
             (unless (flmatrix? M) 
               (raise-type-error 'in-flcolumn "expected matrix, got ~a" M))
             (unless (integer? s) 
               (raise-type-error 'in-flcolumn "expected column number, got ~a" s))
             (unless (and (integer? s) (and (<= 0 s ) (< s n))) 
               (raise-type-error 'in-flcolumn "expected column number, got ~a" s)))
           ([j 0])
           (< j m)
           ([(x) (unsafe-ref a lda j s)])
           #true
           #true
           [(+ j 1)]))]
      [[(i x) (_ M-expr s-expr)]
       #'((x)
          (:do-in
           ([(M s m n a lda) 
             (let ([M1 M-expr])
               (define-param (rd cd a lda) M1)
               (values M1 s-expr rd cd a lda))])
           (begin 
             (unless (flmatrix? M) 
               (raise-type-error 'in-column "expected matrix, got ~a" M))
             (unless (integer? s) 
               (raise-type-error 'in-column "expected col number, got ~a" s))
             (unless (and (integer? s) (and (<= 0 s ) (< s n))) 
               (raise-type-error 'in-column "expected col number, got ~a" s)))
           ([j 0])
           (< j m)
           ([(x) (unsafe-ref a lda j s)]
            [(i) j])
           #true
           #true
           [(+ j 1)]))]
      [[_ clause] (raise-syntax-error 
                   'in-flcolumn "expected (in-flcolumn <flmatrix> <column>)" #'clause #'clause)])))


(define (flmatrix-vandermonde xs n)
  ; One row for each element of xs.
  ; Each row consist of the first 0..n-1 powers of x.
  (define m (length xs))
  (define αs (list->vector xs))
  (define α^j (make-vector m 1.0))
  (for*/flmatrix m n #:column
                 ([j (in-range 0 n)]
                  [i (in-range 0 m)])
                 (define αi^j (vector-ref α^j i))
                 (define αi   (vector-ref αs i ))
                 (vector-set! α^j i (* αi^j αi))
                 αi^j))

;;;
;;; SYNTAX
;;;

(require
 (for-syntax racket/base
             syntax/parse))

(define-syntax (flmatrix: stx)
  (syntax-parse stx 
    [(_ [[x0 xs0 ...] [x xs ...] ...])
     (syntax/loc stx (vectors->flmatrix (vector (vector x0 xs0 ...) (vector x xs ...) ...)))]
    [(_ [xs ... (~and [] r) ys ...])
     (raise-syntax-error 'flmatrix: "given empty row" stx #'r)]
    [(_ (~and [] c))
     (raise-syntax-error 'flmatrix: "given empty matrix" stx #'c)]
    [(_ x)
     (raise-syntax-error 'flmatrix: "expected two-dimensional data" stx)]))

(define-syntax (flrow-matrix stx)
  (syntax-parse stx
    [(_ [x xs ...]) (syntax/loc stx (flmatrix: [[x xs ...]]))]
    [(_ (~and [] r))
     (raise-syntax-error 'flrow-matrix "given empty row" stx #'r)]))

(define-syntax (flcol-matrix stx)
  (syntax-parse stx 
    [(_ [x xs ...])      (syntax/loc stx (flmatrix: [[x] [xs] ...]))]
    [(_ (~and [] c))
     (raise-syntax-error 'flrow-matrix "given empty column" stx #'c)]))


; TODO:
#;(provide 
   ; DONE matrix*
   ; DONE matrix-expt
   ; DONE matrix-ref
   ; DONE matrix-scale
   ; DONE matrix-row-vector?
   ; DONE matrix-column-vector?
   ; DONE matrix/dim     ; construct
   ; DONE matrix-augment ; horizontally
   ; DONE matrix-stack   ; vertically
   ; DONE matrix-block-diagonal
   ; norms
   ; DONE matrix-norm
   ; operators
   ; DONE matrix-transpose
   ; NO-COMPLEX matrix-conjugate
   ; NO-COMPLEX matrix-hermitian
   ; DONE matrix-inverse
   ; row and column
   ; DONE matrix-scale-row
   ; DONE matrix-scale-column
   ; DONE matrix-swap-rows
   ; DONE matrix-swap-columns
   ; DONE matrix-add-scaled-row
   ; DONE ADDED matrix-add-scaled-column
   ; reduction
   ; (DONE) matrix-gauss-eliminate          ; ? use upper in LU ?
   ; DONE matrix-gauss-jordan-eliminate     ; ? LU ? Use matrix-gauss-eliminate
   ; (DONE) matrix-row-echelon-form         ; ? LU ?
   ; DONE matrix-reduced-row-echelon-form    ; ? LU ? Use matrix-gauss-eliminate
   
   ; invariant
   ; DONE matrix-rank  (uses SVD!)
   ; DONE matrix-nullity
   ; DONE matrix-determinant
   ; DONE matrix-trace
   ; spaces
   ;matrix-column+null-space
   ; solvers
   ; DONE matrix-solve
   ; DONE matrix-solve-many
   ; spaces
   matrix-column-space  ; use SVD somehow
   ; column vectors
   ; DONE column        ; construct
   ; DONE unit-column
   ; DONE result-column ; convert to lazy
   ; DONE column-dimension
   ; DONE column-dot
   ; DONE column-norm
   ; DONE column-projection
   ; DONE column-normalize 
   ; DONE scale-column
   ; DONE column+
   ; projection
   ; DONE projection-on-orthogonal-basis
   ; DONE projection-on-orthonormal-basis
   ; DONE projection-on-subspace
   ; DONE gram-schmidt-orthogonal
   ; DONE gram-schmidt-orthonormal
   ; factorization
   ; DONE matrix-lu (renamed to matrix-plu)
   ; DONE matrix-qr
   ; comprehensions
   ; DONE for/matrix:
   ; DONE for*/matrix:
   ; DONE for/matrix-sum:
   ; sequences
   ; DONE in-row
   ; DONE in-column
   ; special matrices
   ; DONE vandermonde-matrix
   )


(define (flmatrix->lists A)
  (map vector->list
       (vector->list 
        (flmatrix->vectors A))))

(define (lists->flmatrix xss)
  (vectors->flmatrix
   (list->vector (map list->vector xss))))


(define (flmatrix-map! A f)
  (define-param (m n a lda) A)
  (for* ([i (in-range m)]
         [j (in-range n)])
    (define aij (unsafe-ref a lda i j))
    (define x   (f aij))
    (define x*  (real->double-flonum x))
    (unsafe-set! a lda i j x*))
  A)


(define (flmatrix-map A f)
  (flmatrix-map! (copy-flmatrix A) f))

(define (flmatrix-make-diagonal v [k 0])
  ; create a square matrix A with diagonal elements from v
  ; if k is given, place the elements on the kth digonal,
  ; k=0 is the main diagonal
  ; k>0 above the main diagonal
  ; k<0 below the main diagonal
  (check-vector  'flmatrix-make-diagonal v)
  (check-integer 'flmatrix-make-diagonal k)
  (define s (+ (vector-length v) (abs k)))
  (define A (make-flmatrix s s))
  (define-param (m n a lda) A)
  (cond
    [(= k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda i i x*))]
    [(> k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda i (+ i k) x*))]
    [(< k 0) (for ([i (in-naturals)] [x (in-vector v)])
               (define x* (real->double-flonum x))
               (unsafe-set! a lda (- i k) i x*))])
  A)

(define (flmatrix-diagonal A [k 0])
  ; extract the k'th diagonal of the matrix A
  (check-flmatrix 'flmatrix-diagonal A)
  (check-integer  'flmatrix-diagonal k)
  (define-param (m n a lda) A)
  (define s (min m n))
  (cond
    [(= k 0) (for/vector ([i (in-range s)])
               (unsafe-ref a lda i i))]
    [(> k 0) (for/vector ([i (in-range (- s k))])
               (unsafe-ref a lda i (+ i k)))]
    [(< k 0) (for/vector ([i (in-range (+ s k))])
               (unsafe-ref a lda i (- i k)))]))

(define (flmatrix-lower-triangle A [k 0])
  ; return triangle with elements on or below the k'th diagonal
  (check-flmatrix 'flmatrix-lower-triangle A)
  (check-integer  'flmatrix-lower-triangle k)
  (define B (copy-flmatrix A))
  (define-param (m n a lda) A)
  (cond
    [(= k 0)  (for* ([i (in-range m)]
                     [j (in-range (+ i 1) n)])
                (flmatrix-set! B i j 0.0))]
    [(> k 0)  (for* ([i (in-range m)]
                     [j (in-range (+ i 1 k) n)])
                (flmatrix-set! B i j 0.0))]
    [(< k 0)  (for* ([i (in-range m)]
                     [j (in-range (max 0 (+ i 1 k)) n)])
                (flmatrix-set! B i j 0.0))])
  B)

(define (flmatrix-circulant-matrix v)
  ; A circulant matrix is a matrix in which each row
  ; is the previous row shifted one to the right.
  (define n (vector-length v))
  (for*/flmatrix n n 
                ([i (in-range n)]
                 [j (in-range n)])
     (vector-ref v (remainder (+ i j) n))))

  
(define (flmatrix-outer-product A B)
  ; accept standard vectors as input
  (define A1 (if (vector? A) (vector->flcolumn A) A))
  (define B1 (if (vector? B) (vector->flrow    B) B))
  ; compute outer product between first column of A and first row of B
  (define-values (am an) (flmatrix-dimensions A1))
  (define-values (bm bn) (flmatrix-dimensions B1))
  (for*/flmatrix am bn ([a (in-flcolumn A1 0)]
                        [b (in-flrow    B1 0)])
                 (* a b)))


; Since the matrix entries of a column are stored contigious,
; we can use cblas_ixamax with incX=1 to find pivots.
(define (flmatrix-find-partial-pivot A i j)
  ; Find the index k of the element a_kj with k>=i
  ; that has the largest absolute value.
  ; I.e. a partial pivot in row j.
  (define-param (m n a lda) A)
  (define ptr (ptr-elm a lda i j)) ; address of the (i,j)th element.
  (define idx (cblas_idamax (- m i) ptr 1))
  (+ i idx))


(define (flmatrix-gauss-elim! A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (define A-original A)
  
  (let loop ([A A] [i 0] [j 0] [without-pivot '()])
    (define-param (m n a lda) A)

    (define (eliminate-row! i pivot l)
      ; eliminate row l using row i which has pivot 
      (define x (unsafe-ref a lda l 0))                   ; first element in row l
      (flmatrix-add-scaled-row! A l (* -1 (/ x pivot)) i) ; scale and subtract
      (unsafe-set! a lda l 0 0.0))                        ; exact 0.0
    
    (define (eliminate-rows-below! i pivot) (for ([l (in-range (+ i 1) m)]) (eliminate-row! i pivot l)))
    (define (eliminate-rows-above! i pivot) (for ([l (in-range 0 i)])       (eliminate-row! i pivot l)))

    (define (A-without-first-column) (shared-submatrix! A 0 1 m (- n 1)))
    
    (cond
      [(= n 0) (values A-original (reverse without-pivot))]
      ;; None of the rest of the columns can have pivots
      [(= i m) (values A-original (append (reverse without-pivot) (range j (+ j n))))]
      [else    (define p  (case pivoting
                            [(partial) (flmatrix-find-partial-pivot A i 0)]
                            #;[(first)   (flmatrix-find-first-pivot   A 0 0)]
                            [else (error 'flmatrix-gauss-elim! "unknown pivoting type")]))
               (define pivot (flmatrix-ref A p 0))
               (cond                 
                 [(<= pivot epsilon) ;; no pivot
                  (loop (A-without-first-column) i (+ j 1) (cons j without-pivot))]
                 [else               ;; pivot found
                  (flmatrix-swap-rows! A i p)
                  (eliminate-rows-below! i pivot)
                  (when jordan?        (eliminate-rows-above! i pivot))
                  (when unitize-pivot? (flmatrix-scale-row! A i (/ 1. pivot)))
                  (loop (A-without-first-column) (+ i 1) (+ j 1) without-pivot)])])))

(define (flmatrix-gauss-elim A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (flmatrix-gauss-elim! (copy-flmatrix A) jordan? unitize-pivot? pivoting))

(define (matrix-row-echelon! A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (flmatrix-gauss-elim! A jordan? unitize-pivot? pivoting)
  A)

(define (matrix-row-echelon A [jordan? #f] [unitize-pivot? #f] [pivoting 'partial])
  (matrix-row-echelon! (copy-flmatrix A) jordan? unitize-pivot? pivoting))




;;;
;;; TEST
;;;


(module+ test
  (require rackunit)
  
  (define (flcheck-equal? a b)
    (< (abs (- b a)) 0.00001))
  

  (with-check-info
   (['test-case "flmatrix/dim"])
   (check-equal? (flmatrix->vector (flmatrix/dim 2 2 1 2 3 4))
                 #(1. 2. 3. 4.))
   (check-equal? (flmatrix->vector (vector->flmatrix 2 2 #(1 2 3 4)))
                 #(1. 2. 3. 4.))
   (check-equal? (flmatrix->vector (flmatrix/dim 2 2 1 2 3 4))
                 #(1. 2. 3. 4.))
   (check-equal? (flmatrix->vectors (flmatrix/dim 2 2 1 2 3 4))
                 #(#[1. 2.] #[3. 4.]))

   (let ()
     (define A  (flmatrix/dim 2 2 1 2 3 4))
     (define B  (flmatrix/dim 2 2 5 6 7 8))
     (define AB (flmatrix/dim 2 2 19 22 43 50))
     (check-equal? (flmatrix->vectors (flmatrix* A B))
                   (flmatrix->vectors AB)))

   (let ()
     (define C  (flmatrix/dim 2 2 1 2 3 4))
     (define D  (flmatrix/dim 2 3 5 6 7 8 9 10))
     (define CD (flmatrix/dim 2 3 21 24 27 47 54 61))
     (check-equal? (flmatrix->vectors (flmatrix* C D))
                   (flmatrix->vectors CD)))

   (check-equal? (flmatrix->vectors 
                  (flmatrix* (flmatrix/dim 2 3 0 0 1 0 0 0)
                             (flmatrix/dim 3 2 1 2 3 4 5 6)))
                 (flmatrix->vectors 
                  (flmatrix/dim 2 2 5 6 0 0))))


  (with-check-info
   (['test-group "matrix-constructors.rkt"])
   (with-check-info
    (['test-case 'flmatrix-identity])
          
    (check-equal? (flmatrix->lists (flmatrix-identity 1)) '[[1.]])
    (check-equal? (flmatrix->lists (flmatrix-identity 2)) '[[1. 0.] [0. 1.]])
    (check-equal? (flmatrix->lists (flmatrix-identity 3)) '[[1. 0. 0.] [0. 1. 0.] [0. 0. 1.]]) 
    (check-equal? (flmatrix->lists (flmatrix-identity 1)) '[[1.]])
    (check-equal? (flmatrix->lists (flmatrix-identity 2)) '[[1. 0.] [0. 1.]])
    (check-equal? (flmatrix->lists (flmatrix-identity 3)) '[[1. 0. 0.] [0. 1. 0.] [0. 0. 1.]]))
   (with-check-info
    (['test-case 'const-matrix])          
    (check-equal? (flmatrix->lists (make-flmatrix 2 3 0.)) '((0. 0. 0.) (0. 0. 0.))))
   (with-check-info
    (['test-case 'matrix->list])
    (check-equal? (flmatrix->lists (lists->flmatrix '((1. 2.) (3. 4.)))) '((1. 2.) (3. 4.))))
   (with-check-info
    (['test-case 'matrix->vectors])
    (check-equal? (flmatrix->vectors (vectors->flmatrix '#(#(1. 2.) #(3. 4.)))) '#(#(1. 2.) #(3. 4.))))
   (with-check-info
    (['test-case 'matrix-row])
    (check-equal? (flmatrix-row (flmatrix-identity 3) 0) (list->flmatrix '[[1 0 0]]))
    (check-equal? (flmatrix-row (flmatrix-identity 3) 1) (list->flmatrix '[[0 1 0]]))
    (check-equal? (flmatrix-row (flmatrix-identity 3) 2) (list->flmatrix '[[0 0 1]])))
   (with-check-info
    (['test-case 'matrix-col])
    (check-equal? (flmatrix-column (flmatrix-identity 3) 0) (list->flmatrix '[[1] [0] [0]]))
    (check-equal? (flmatrix-column (flmatrix-identity 3) 1) (list->flmatrix '[[0] [1] [0]]))
    (check-equal? (flmatrix-column (flmatrix-identity 3) 2) (list->flmatrix '[[0] [0] [1]])))
   (with-check-info
    (['test-case 'flsubmatrix])
    (check-equal? (flsubmatrix (flmatrix-identity 3) 1 2 0 0)
                  (list->flmatrix '[[1 0]]))
    (check-equal? (flsubmatrix (flmatrix-identity 3) 2 3 0 0)
                  (list->flmatrix '[[1 0 0] [0 1 0]]))))

  (with-check-info
    (['test-group "product"])
    (check-equal?
     (flmatrix*vector (list->flmatrix '[[1 2 3] [4 5 6]])
                      (list->flmatrix '[[10] [11] [12]]))
     (list->flmatrix '[[68] [167]]))
    (check-equal?
     (flmatrix*vector (list->flmatrix '[[1 4] [2 5] [3 6]])
                      (list->flmatrix '[[10] [11] [12]])
                      #f 1. 1. #t)
     (list->flmatrix '[[68] [167]])))

  (with-check-info
   (['test-group "flmatrix-pointwise.rkt"])
   (let ()
     (define A   (list->flmatrix '[[1 2] [3 4]]))
     (define ~A  (list->flmatrix '[[-1 -2] [-3 -4]]))
     (define B   (list->flmatrix '[[5 6] [7 8]]))
     (define A+B (list->flmatrix '[[6 8] [10 12]]))
     (define A-B (list->flmatrix '[[-4 -4] [-4 -4]]))         
     (with-check-info
      (['test-case 'flmatrix+])
      (check-equal? (flmatrix+ A B) A+B))
     (with-check-info
      (['test-case 'flmatrix-])
      (check-equal? (flmatrix- A B) A-B)
      (check-equal? (flmatrix- A)   ~A))))
  
  (with-check-info
   (['test-group "flmatrix-expt.rkt"])
   (define A (list->flmatrix '[[1 2] [3 4]]))
   (with-check-info
    (['test-case 'flmatrix-expt])
    (check-equal? (flmatrix-expt A 0) (flmatrix-identity 2))
    (check-equal? (flmatrix-expt A 1) A)
    (check-equal? (flmatrix-expt A 2) (list->flmatrix '[[7 10] [15 22]]))
    (check-equal? (flmatrix-expt A 3) (list->flmatrix '[[37 54] [81 118]]))
    (check-equal? (flmatrix-expt A 8) (list->flmatrix '[[165751 241570] [362355 528106]]))))

  (with-check-info
   (['test-group "flmatrix-operations.rkt"])
   (with-check-info
    (['test-case 'vandermonde-flmatrix])
    (check-equal? (flmatrix-vandermonde '(1 2 3) 5)
                  (list->flmatrix '[[1 1 1 1 1] [1 2 4 8 16] [1 3 9 27 81]])))
   (with-check-info
    (['test-case 'in-column])
    (check-equal? (for/list ([x (in-flcolumn (flmatrix/dim 2 2  1 2 3 4) 0)]) x)
                  '(1. 3.))
    (check-equal? (for/list ([x (in-flcolumn (flmatrix/dim 2 2  1 2 3 4) 1)]) x)
                  '(2. 4.))
    (check-equal? (for/list ([x (in-flcolumn (flcolumn 5 2 3))]) x)
                  '(5. 2. 3.)))
   (with-check-info
    (['test-case 'in-row])
    (check-equal? (for/list ([x (in-flrow (flmatrix/dim 2 2  1 2 3 4) 0)]) x)
                  '(1. 2.))
    (check-equal? (for/list ([x (in-flrow (flmatrix/dim 2 2  1 2 3 4) 1)]) x)
                  '(3. 4.)))
   (with-check-info
    (['test-case 'for/flmatrix:])
    (check-equal? (for/flmatrix 2 4 ([i (in-naturals)]) i)
                  (flmatrix/dim 2 4 
                                0 1 2 3
                                4 5 6 7))
    (check-equal? (for/flmatrix 2 4 #:column ([i (in-naturals)]) i)
                  (flmatrix/dim 2 4    
                                0 2 4 6
                                1 3 5 7))
    (check-equal? (for/flmatrix 3 3 ([i (in-range 10 100)]) i)
                  (flmatrix/dim 3 3 10 11 12 13 14 15 16 17 18)))
   (with-check-info
    (['test-case 'for*/flmatrix:])
    (check-equal? (for*/flmatrix 3 3 ([i (in-range 3)] [j (in-range 3)]) (+ (* i 10) j))
                  (flmatrix/dim 3 3 0 1 2 10 11 12 20 21 22)))    
   (with-check-info
    (['test-case 'flmatrix-block-diagonal])
    (check-equal? (flmatrix-block-diagonal (flmatrix/dim 2 2 1 2 3 4) (flmatrix/dim 1 3 5 6 7))
                  (list->flmatrix '[[1 2 0 0 0] [3 4 0 0 0] [0 0 5 6 7]])))
   (with-check-info
    (['test-case 'flmatrix-augment])
    (check-equal? (flmatrix-augment (flcolumn 1 2 3) (flcolumn 4 5 6) (flcolumn 7 8 9))
                  (flmatrix/dim 3 3  1 4 7  2 5 8  3 6 9)))
   (with-check-info
    (['test-case 'flmatrix-stack])
    (check-equal? (flmatrix-stack (flcolumn 1 2 3) (flcolumn 4 5 6) (flcolumn 7 8 9))
                  (flcolumn 1 2 3 4 5 6 7 8 9)))
   (with-check-info
    (['test-case 'column-dimension])
    (= (flcolumn-size #(1 2 3)) 3)
    (= (flcolumn-size (vector->flmatrix 1 2 #(1 2))) 1))
   (let ([flmatrix: vector->flmatrix])
     (with-check-info
      (['test-case 'column-dot])
      (= (flcolumn-dot (flcolumn 1 2)   (flcolumn 1 2)) 5)
      (= (flcolumn-dot (flcolumn 1 2)   (flcolumn 3 4)) 11)
      (= (flcolumn-dot (flcolumn 3 4)   (flcolumn 3 4)) 25)
      (= (flcolumn-dot (flcolumn 1 2 3) (flcolumn 4 5 6))
         (+ (* 1 4) (* 2 5) (* 3 6)))))
   (with-check-info
    (['test-case 'flmatrix-trace])
    (check-equal? (flmatrix-trace (vector->flmatrix 2 2 #(1 2 3 4))) 5.))
   (let ([flmatrix: vector->flmatrix])
     (with-check-info
      (['test-case 'column-norm])
      (= (flcolumn-norm (flcolumn 2 4)) (sqrt 20))))
   (with-check-info
    (['test-case 'column-projection])
    (check-equal? (flcolumn-projection #(1 2 3) #(4 5 6)) (flcolumn 128/77 160/77 192/77))
    (check-equal? (flcolumn-projection (flcolumn 1 2 3) (flcolumn 2 4 3))
                  (flmatrix-scale 19/29 (flcolumn 2 4 3))))
   (with-check-info
    (['test-case 'projection-on-orthogonal-basis])
    (check-equal? (flprojection-on-orthogonal-basis #(3 -2 2) (list #(-1 0 2) #( 2 5 1)))
                  (flcolumn -1/3 -1/3 1/3))
    (check-equal? (flprojection-on-orthogonal-basis 
                   (flcolumn 3 -2 2) (list #(-1 0 2) (flcolumn 2 5 1)))
                  (flcolumn -1/3 -1/3 1/3)))
   (with-check-info
    (['test-case 'projection-on-orthonormal-basis])
    (check-equal? (flprojection-on-orthonormal-basis 
                   #(1 2 3 4) 
                   (list (flmatrix-scale 1/2 (flcolumn  1  1  1 1))
                         (flmatrix-scale 1/2 (flcolumn -1  1 -1 1))
                         (flmatrix-scale 1/2 (flcolumn  1 -1 -1 1))))
                  (flcolumn 2 3 2 3)))
   (with-check-info
    (['test-case 'flgram-schmidt-orthogonal])
    (check-equal? (flgram-schmidt-orthogonal (list #(3 1) #(2 2)))
                  (list (flcolumn 3 1) (flcolumn -2/5 6/5))))
   (with-check-info
    (['test-case 'flvector-normalize])
    (check-equal? (flcolumn-normalize #(3 4)) 
                  (flcolumn 3/5 4/5)))
   (with-check-info
    (['test-case 'flgram-schmidt-orthonormal])
    (check-equal? (flgram-schmidt-orthonormal '(#(3 1) #(2 2)))
                  (list (flcolumn-normalize #(3 1))
                        (flcolumn-normalize #(-2/5 6/5)))))
    
   (with-check-info
    (['test-case 'projection-on-subspace])
    (check-equal? (flprojection-on-subspace #(1 2 3) '(#(2 4 3)))
                  (flmatrix-scale 19/29 (flcolumn 2 4 3))))
   (with-check-info
    (['test-case 'unit-vector])
    (check-equal? (flcolumn-unit 4 1) (flcolumn 0 1 0 0)))
    (with-check-info (['test-case 'flmatrix-qr])
      (let*-values ([(A) (flmatrix/dim 3 2  1 1 0 1 1 1)]
                    [(Q R) (flmatrix-qr A)])
        (check-true
         (flmatrix= (flmatrix* Q R)
                    A
                    epsilon))))
   (with-check-info
    (['test-case 'flmatrix-solve])
    (let* ([M (list->flmatrix '[[1 5] [2 3]])] 
           [b (list->flmatrix '[[5] [5]])])
      (check-equal? (flmatrix* M (flmatrix-solve M b)) b)))
   (with-check-info
    (['test-case 'flmatrix-inverse])
    (check-equal? (let ([M (list->flmatrix '[[1 2] [3 4]])]) (flmatrix* M (flmatrix-inverse M)))
                  (flmatrix-identity 2))
    (check-equal? (let ([M (list->flmatrix '[[1 2] [3 4]])]) (flmatrix* (flmatrix-inverse M) M))
                  (flmatrix-identity 2)))
   (with-check-info
    (['test-case 'flmatrix-determinant])
    (check-equal? (flmatrix-determinant (list->flmatrix '[[3]])) 3.)
    (check-equal? (flmatrix-determinant (list->flmatrix '[[1 2] [3 4]])) (- (* 1. 4.) (* 2. 3.)))
    (flcheck-equal? (flmatrix-determinant (list->flmatrix '[[1 2 3] [4  5 6] [7 8 9]])) 0.)
    (flcheck-equal? (flmatrix-determinant (list->flmatrix '[[1 2 3] [4 -5 6] [7 8 9]])) 120.)
    (flcheck-equal? (flmatrix-determinant 
                     (list->flmatrix '[[1 2 3 4] [-5 6 7 8] [9 10 -11 12] [13 14 15 16]])) 5280.))
   (with-check-info
    (['test-case 'flmatrix-scale])
    (check-equal? (flmatrix-scale 2 (list->flmatrix '[[1 2] [3 4]]))
                  (list->flmatrix '[[2 4] [6 8]])))
   (with-check-info
    (['test-case 'flmatrix-transpose])
    (check-equal? (flmatrix-transpose (list->flmatrix '[[1 2] [3 4]]))
                  (list->flmatrix '[[1 3] [2 4]])))
   ; TODO: Just use U from LU factorization
   #;(let ()
       (: gauss-eliminate : (flmatrix Number) Boolean Boolean -> (flmatrix Number))
       (define (gauss-eliminate M u? p?)
         (let-values ([(M wp) (flmatrix-gauss-eliminate M u? p?)])
           M))
       (with-check-info
        (['test-case 'flmatrix-gauss-eliminate])
        (check-equal? (let ([M (list->flmatrix '[[1 2] [3 4]])])
                        (gauss-eliminate M #f #f))
                      (list->flmatrix '[[1 2] [0 -2]]))
        (check-equal? (let ([M (list->flmatrixixix  '[[2 4] [3 4]])])
                        (gauss-eliminate M #t #f))
                      (list->flmatrixixixix '[[1 2] [0 1]]))
        (check-equal? (let ([M (list->flmatrixix  '[[2. 4.] [3. 4.]])])
                        (gauss-eliminate M #t #t))
                      (list->flmatrixix '[[1. 1.3333333333333333] [0. 1.]]))
        (check-equal? (let ([M (list->flmatrix  '[[1 4] [2 4]])])
                        (gauss-eliminate M #t #t))
                      (list->flmatrix '[[1 2] [0 1]]))
        (check-equal? (let ([M (list->flmatrix  '[[1 2] [2 4]])])
                        (gauss-eliminate M #f #t))
                      (list->flmatrix '[[2 4] [0 0]]))))
   (with-check-info
    (['test-case 'flmatrix-scale-row])
    (check-equal? (flmatrix-scale-row (flmatrix-identity 3) 0 2)
                  (lists->flmatrix '[[2 0 0] [0 1 0] [0 0 1]])))
   (with-check-info
    (['test-case 'flmatrix-swap-rows])
    (check-equal? (flmatrix-swap-rows (lists->flmatrix '[[1 2 3] [4 5 6] [7 8 9]]) 0 1)
                  (lists->flmatrix '[[4 5 6] [1 2 3] [7 8 9]])))
   (with-check-info
    (['test-case 'flmatrix-add-scaled-row])
    (check-equal? (flmatrix-add-scaled-row (lists->flmatrix '[[1 2 3] [4 5 6] [7 8 9]]) 0 2 1)
                  (lists->flmatrix '[[9 12 15] [4 5 6] [7 8 9]])))
   (let ()
     (define M (lists->flmatrix '[[1  1  0  3]
                                  [2  1 -1  1]
                                  [3 -1 -1  2]
                                  [-1  2  3 -1]]))
     (define-values (P L U) (flmatrix-plu M))
     (with-check-info
      (['test-case 'flmatrix-plu])
      (check-equal? (flmatrix* P (flmatrix* L U)) M)))
   (with-check-info
    (['test-case 'flmatrix-rank])
    (check-equal? (flmatrix-rank (list->flmatrix '[[0 0] [0 0]])) 0)
    (check-equal? (flmatrix-rank (list->flmatrix '[[1 0] [0 0]])) 1)
    (check-equal? (flmatrix-rank (list->flmatrix '[[1 0] [0 3]])) 2)
    (check-equal? (flmatrix-rank (list->flmatrix '[[1 2] [2 4]])) 1)
    (check-equal? (flmatrix-rank (list->flmatrix '[[1 2] [3 4]])) 2))
   (with-check-info
    (['test-case 'flmatrix-nullity])
    (check-equal? (flmatrix-nullity (list->flmatrix '[[0 0] [0 0]])) 2)
    (check-equal? (flmatrix-nullity (list->flmatrix '[[1 0] [0 0]])) 1)
    (check-equal? (flmatrix-nullity (list->flmatrix '[[1 0] [0 3]])) 0)
    (check-equal? (flmatrix-nullity (list->flmatrix '[[1 2] [2 4]])) 1)
    (check-equal? (flmatrix-nullity (list->flmatrix '[[1 2] [3 4]])) 0))
   ; Not implemented yet...
   #;(let ()
       (define-values (c1 n1) 
         (flmatrix-column+null-space (list->flmatrix '[[0 0] [0 0]])))
       (define-values (c2 n2) 
         (flmatrix-column+null-space (list->flmatrix '[[1 2] [2 4]])))
       (define-values (c3 n3) 
         (flmatrix-column+null-space (list->flmatrix '[[1 2] [2 5]])))
       (with-check-info
        (['test-case 'flmatrix-column+null-space])
        (check-equal? c1 '())
        (check-equal? n1 (list (list->flmatrix '[[0] [0]])
                               (list->flmatrix '[[0] [0]])))
        (check-equal? c2 (list (list->flmatrix '[[1] [2]])))
        ;(check-equal? n2 '([0 0]))
        (check-equal? c3 (list (list->flmatrix '[[1] [2]])
                               (list->flmatrix '[[2] [5]])))
        (check-equal? n3 '()))))
  
  (with-check-info
   (['test-group "matrix-multiply.rkt"])
   (with-check-info
    (['test-case 'flmatrix*])
    (let ()
      (define-values (A B AB) (values '[[1 2] [3 4]] '[[5 6] [7 8]] '[[19 22] [43 50]]))
      (check-equal? (flmatrix* (list->flmatrix A) (list->flmatrix B)) (list->flmatrix AB)))
    (let () 
      (define-values (A B AB) (values '[[1 2] [3 4]] '[[5 6 7] [8 9 10]] '[[21 24 27] [47 54 61]]))
      (check-equal? (flmatrix* (list->flmatrix A) (list->flmatrix B)) (list->flmatrix AB))))))

(define (build-flmatrix m n f)
  (for*/flmatrix m n 
                 ([i (in-range m)]
                  [j (in-range n)])
                 (f i j)))

