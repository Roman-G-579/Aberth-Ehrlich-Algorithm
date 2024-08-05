#lang racket

(define (round-complex complex-nums)
  (map (λ(num)
         (let ((real (real-part num))
               (imag (imag-part num)))
           (+ (/ (round (* real 1000)) 1000)
              (* (/ (round (* imag 1000)) 1000) 0+i))))
       complex-nums))

; Reads the file line-by-line, trims each line and converts the string to a number,
; that is added to a list
; Once end-of-file is reached, the list is returned
(define (read-coefficients file)
  (let ((line (read-line file 'any)))
    (if (eof-object? line)
        '()
         (cons (string->number (string-trim line)) (read-coefficients file)))))


(define (is-close-to-zero value epsilon)
  (and (< (abs (real-part value)) epsilon)
       (< (abs (imag-part value)) epsilon)))

(define pi 3.141592653589793)

 (define poly-val (λ(p x)
         (if (null? p)
             0
             (+ (car p) (* x (poly-val (cdr p) x))))))

(define (frac-val p q x)
  (if (<= (magnitude x) 1)
      (/ (poly-val (reverse p) x) (poly-val (reverse q) x))
      (* x (/ (poly-val p (/ 1 x))
              (poly-val q (/ 1 x))))))

(define (polynomial-derivative-coefficients coefficients)
  (let ((n (length coefficients))) ; n size is equal to the amount of coefficients in the polynomial
    (if (= n 1)
        '(0)
        (let loop ((i 0) ; loop variable is initialized to 0
                   (result '()))
          (if (= i (- n 1))
              (reverse result) ; If all elements have been processed, reverse the list

              ; Recursively calls loop with the next index
              ; and the current result list
              (loop (+ i 1)
                    (cons (* (list-ref coefficients i) (- n i 1)) result))))))) ; Calculates the coeff of the i-th term

(define (get-approximations coefficients)
  (let* ((deg (- (length coefficients) 1))
         (pi-val pi)
         (p-0 (list-ref coefficients deg))
         (p-n (car coefficients))
         (r (magnitude (expt (/ p-0 p-n) (/ 1 deg))))
         (angles (for/list ((x (in-range 0 deg))) (* x (/ (* 2 pi-val) deg)))))
    (map (λ(angle) (+ (* r (cos angle)) (* r (sin angle) 0+i))) angles)))


(define (calculate-sigma approximations k)
  (let ((zk (list-ref approximations k)))
    (foldl (λ(zj sigma)
             (if (not (= zj zk))
                 (+ sigma (/ 1 (- zk zj)))
                 sigma))
           0.0
           approximations)))

; Check if all given approximations count as converged for the given epsilon
; by plugging each approximation into the polynomial and evaluating the result's distance from the epsilon
; If all converge, returns true, otherwise returns false
(define (check-convergence approximations epsilon)
  (let* ((rev-coeffs (reverse coefficients)) ; Coeffs are reversed in advance before being sent to the poly-val function
         (poly-evals (map (λ(val) (poly-val rev-coeffs val)) approximations))) ; All polynomial values are precomputed:
                                                                                ; poly-val is applied on every element
                                                                                ; of the approximations list using
                                                                                ; the lambda function
    ; Every element in poly-eval is sent as a parameter to is-close-to-zero,
    ; along with the epsilon.
    ; If any element isn't precise enough, returns false
    ; If all elements are precise enough, returns true
    (foldl (λ(eval roots-converge)
                (if (not (is-close-to-zero eval epsilon))
                    #f
                    roots-converge))
           #t ; The initial value of roots-converge
           poly-evals))) ; List of polynomial evaluations

(define (aberth-method coefficients epsilon)
  (let ((approximations (get-approximations coefficients))
         (derivative-coefficients (polynomial-derivative-coefficients coefficients))
         (max-iterations 20))

        (define (iterate approximations iteration)
          (if (>= iteration max-iterations)
              approximations ; If max-iterations has been reached, the approximations list is returned
              (let ((offsets
                     (for/list ((zk approximations) (k (in-naturals)))
                       (let* ((frac (frac-val coefficients derivative-coefficients zk))
                              (sigma (calculate-sigma approximations k)))
                         (/ frac (- 1 (* frac sigma)))))))

             (let ((new-approximations (map - approximations offsets))) ;

             (if (check-convergence new-approximations epsilon)
                new-approximations ; If all new approximations are sufficiently precise, returns them as a list

        (iterate new-approximations (+ iteration 1)))))))

    (iterate approximations 0)) ; Initial call to iterate
)


; Opens an input file (calls open-input-file), and calls read-coefficients using it as a parameter
(define coefficients (call-with-input-file "poly_coeff(997).txt" read-coefficients))

(define result (time (aberth-method coefficients 0.0001)))


;(printf "Coefficients: ~a\n" derivative-coefficients)
;(printf "Roots: ~a\n" result)
;(printf "Rounded Roots: ~a\n" (round-complex result))
;(printf "len ~a\n" (length result))

; Define the function to print each complex number on a new line
(define (print-complex-list complex-list)
  (for-each (lambda (z) (printf "~a\n" z)) complex-list))

(print-complex-list result)