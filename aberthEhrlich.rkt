#lang racket

(define Factorial_ (λ(n acc)
    (if (= n 0)
    acc
    (Factorial_ (- n 1) (* n acc)))))

(define Factorial (λ(n)(Factorial_ n 1)))

(display (Factorial 6) )
