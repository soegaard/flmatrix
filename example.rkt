#lang racket
(require "flmatrix.rkt")
(define dim 1000)
(define big1 (build-flmatrix dim dim (lambda (i j) (random))))
(define big2 (build-flmatrix dim dim (lambda (i j) (random))))
(define res (time (flmatrix* big1 big2)))