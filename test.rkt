#lang racket
(define Poly_val (λ(p x)
         (if (null? p)
             0
             (+ (car p) (* x (Poly_val (cdr p) x)))
)
))
(define p '(1 2 3 4))

(printf "~a\n" (Poly_val p 2))