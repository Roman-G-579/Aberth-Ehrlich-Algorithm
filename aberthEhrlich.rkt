#lang racket

(define pi 3.141592653589793)


; Reads the file line-by-line, trims each line and converts the string to a number,
; that is added to a list
; Once end-of-file is reached, the list is returned
(define (read-coefficients file)
  (let ((line (read-line file)))
    (if (eof-object? line)
        '()
         (cons (string->number (string-trim line)) (read-coefficients file)))))

; Checks if a complex number is close to zero within a given epsilon
(define (is-close-to-zero value epsilon)
  (and (< (abs (real-part value)) epsilon)
       (< (abs (imag-part value)) epsilon)))

; Evaluates a polynomial at a given value x
; a_0+a_1x+...+a_nx^n = a0 + x(a_1+x(a_2+x(a_3+...+x(a_(n-1)+xa_n)...)))
(define poly-val (λ(p x)
         (if (null? p)
             0 ; Base case: if the list of coefficients is empty, return 0
             ; Otherwise, take the first coefficient, multiply it by x and add the result of evaluating the rest of the polynomial
             (+ (car p) (* x (poly-val (cdr p) x))))))

; Evaluates a rational function p(x)/q(x) at x
(define (frac-val p q x)
  (if (<= (magnitude x) 1)
      (/ (poly-val (reverse p) x) (poly-val (reverse q) x))
      (* x (/ (poly-val p (/ 1 x))
              (poly-val q (/ 1 x))))))


; Takes a list of coefficients, and returns a list of the derivative's coefficients
(define (polynomial-derivative-coefficients coefficients)
  (let ((n (length coefficients))) ; n size is equal to the amount of coefficients in the polynomial
    (if (= n 1)
        '(0)
        (let loop ((i 0) (result '())) ; loop function with i's value initialized to 0
          (if (= i (- n 1))
              (reverse result) ; If all elements have been processed, reverse the list

              ; Recursively calls loop with the next index
              ; and the current result list, with
              (loop (+ i 1) (cons (* (list-ref coefficients i) (- n i 1)) result))))))) ; Calculates the coeff of the i-th term
              ; [4x^2,2x,1] -> 4x2+2x1 -> [8,2]

; Generates initial approximations for the polynomial roots
(define (get-approximations coefficients)
  (let* ((deg (- (length coefficients) 1))
         (p-0 (list-ref coefficients deg))
         (p-n (car coefficients))
         (r (magnitude (expt (/ p-0 p-n) (/ 1 deg)))) ; r = (|p0/pn|)^(1/deg)

         ;Distributes angles evenly, the angle amount is equal to the degree of the polynomial
         (angles (for/list ((x (in-range deg))) (* x (/ (* 2 pi) deg))))) ; Computes the angles for the approximations
    (map (λ(angle) (+ (* r (cos angle)) (* r (sin angle) 0+i))) angles))) ; Converts the polar coordinates to complex numbers


; Calculates the distance from approximation k to every other approximation in the list
(define (calculate-sigma approximations zk)
  (let loop ((lst approximations) (sigma 0.0))
      (cond
        ((null? lst) sigma) ; If the list is empty, return the accumulated sigma
        ((not (= (car lst) zk)) ; If the current element is not equal to zk
         (loop (cdr lst) (+ sigma (/ 1 (- zk (car lst)))))) ; Update sigma and continue with the rest of the list
        (else (loop (cdr lst) sigma))))) ; If the current element equals zk, skip it and continue


; Check if all given approximations count as converged for the given epsilon
; by plugging each approximation into the polynomial and evaluating the result's distance from the epsilon
; If all converge, returns true, otherwise returns false
(define (check-convergence approximations epsilon)
  (let* ((rev-coeffs (reverse coefficients)) ; Coeffs are reversed in advance before being sent to the poly-val function
         (poly-evals (map (λ(val) (poly-val rev-coeffs val)) approximations))) ; All polynomial values are precomputed:
                                                                                ; poly-val is applied on every element
                                                                                ; of the approximations list using
                                                                                ; the lambda function

    (let loop ((lst poly-evals))
        (cond
          ((null? lst) #t) ; If the list is empty, all elements were close to zero
          ((not (is-close-to-zero (car lst) epsilon)) #f) ; If an element is not close to zero, return #f
          (else (loop (cdr lst))))))) ; Otherwise, continue with the rest of the list

; Implements the Aberth method for finding polynomial roots
(define (aberth-method coefficients epsilon)
  (let ((approximations (get-approximations coefficients))
         (derivative-coefficients (polynomial-derivative-coefficients coefficients))
         (max-iterations 20)) ; Constant value for max iteration limit

        (define (iterate approximations iteration)
          (if (>= iteration max-iterations)
              approximations ; If max-iterations has been reached, the approximations list is returned
              (let ((offsets
                     (for/list ((zk approximations))
                       (let ((frac (frac-val coefficients derivative-coefficients zk))
                              (sigma (calculate-sigma approximations zk)))
                         (/ frac (- 1 (* frac sigma))))))) ; offset = frac / ( 1 - (frac*sigma) )

             (let ((new-approximations (map - approximations offsets))) ;

             (if (check-convergence new-approximations epsilon)
                new-approximations ; If all new approximations are sufficiently precise, returns them as a list

        (iterate new-approximations (+ iteration 1)))))))

    (iterate approximations 0)) ; Initial call to iterate
)


; Opens an input file (calls open-input-file), and calls read-coefficients using it as a parameter
(define coefficients (call-with-input-file "poly_coeff(997).txt" read-coefficients))


(define result (time (aberth-method coefficients 0.0001)))


; Prints each complex number in a new line
(define (print-complex-list complex-list)
  (for-each (lambda (z) (printf "~a\n" z)) complex-list))

(print-complex-list result)





;(define test-val (poly-val (reverse coefficients) 0.7450737533631855-0.6700384058892166i))
;(print test-val)