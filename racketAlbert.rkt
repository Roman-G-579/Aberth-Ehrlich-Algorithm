#lang racket

;; This function rounds the real and imaginary parts of a list of complex numbers
(define (round-complex complex-nums)
  (map (lambda (num)                      ; map applies the following lambda function to each element in complex-nums
         (let ((real (real-part num))     ; Extract the real part of the complex number
               (imag (imag-part num)))    ; Extract the imaginary part of the complex number
           ;; Round both real and imaginary parts to three decimal places
           (+ (/ (round (* real 1000)) 1000)
              (* (/ (round (* imag 1000)) 1000) 0+i)))) ; Combine the rounded parts into a new complex number
       complex-nums))                     ; Apply the lambda function to all elements in complex-nums

;; This function reads coefficients from a file
(define (read-coefficients filename)
  (with-input-from-file filename          ; Open the file for reading
    (lambda ()                            ; Define a lambda function to process the file
      (let loop ((lines '()))             ; Start a recursive loop with an empty list of lines
        (let ((line (read-line)))         ; Read a line from the file
          (if (eof-object? line)          ; Check if end of file is reached
              (reverse lines)             ; If so, return the list of lines in reverse order
              (let ((num (string->number (string-trim line)))) ; Convert the line to a number
                (if (number? num)         ; Check if the conversion was successful
                    (loop (cons num lines)) ; Add the number to the list and continue looping
                    (begin
                      (printf "Error: Invalid number format in line: ~a\n" line) ; Print an error message
                      (loop lines)))))))))) ; Continue looping without adding the invalid line

;; This function checks if a complex number is close to zero within a given epsilon
(define (is-close-to-zero? value epsilon)
  (and (< (abs (real-part value)) epsilon) ; Check if the real part is close to zero
       (< (abs (imag-part value)) epsilon))) ; Check if the imaginary part is close to zero

;; Define the value of pi
(define pi 3.141592653589793)

;; This function evaluates a polynomial at a given value x
(define poly-val (λ(p x)
         (if (null? p)                    ; Base case: if the list of coefficients is empty
             0                            ; Return 0
             (+ (car p)                   ; Otherwise, take the first coefficient
                (* x (poly-val (cdr p) x)))))) ; Multiply it by x and add the result of evaluating the rest of the polynomial

;; This function evaluates a rational function p(x)/q(x) at x
(define (frac-val p q x)
  (if (<= (magnitude x) 1)                ; If the magnitude of x is less than or equal to 1
      (/ (poly-val (reverse p) x) (poly-val (reverse q) x)) ; Evaluate p and q at x and return their ratio
      (* x (/ (poly-val p (/ 1 x))        ; Otherwise, evaluate p and q at 1/x
              (poly-val q (/ 1 x))))))    ; Multiply the result by x

;; This function computes the derivative of a polynomial's coefficients
(define (polynomial-derivative-coefficients coefficients)
  (let* ((n (length coefficients))        ; Get the number of coefficients
         (derivative-coeffs '()))         ; Initialize an empty list for the derivative coefficients
    (if (= n 1)                           ; If there is only one coefficient
        '(0)                              ; The derivative is zero
        (for/list ((i (in-range 0 (- n 1)))) ; Otherwise, loop over the coefficients
          (* (- n 1 i)                    ; Multiply each coefficient by its position in the list
             (list-ref coefficients i)))))) ; And store the result in the list of derivative coefficients

;; This function generates initial approximations for the polynomial roots
(define (get-approximations coefficients)
  (let* ((deg (- (length coefficients) 1)) ; Get the degree of the polynomial
         (pi-val pi)                      ; Use the value of pi
         (p-0 (list-ref coefficients deg)) ; Get the leading coefficient
         (p-n (car coefficients))         ; Get the constant term
         (r (magnitude (expt (/ p-0 p-n) (/ 1 deg)))) ; Compute the radius for the approximations
         (angles (for/list ((x (in-range 0 deg))) (* x (/ (* 2 pi-val) deg))))) ; Compute the angles for the approximations
    (map (lambda (angle) (+ (* r (cos angle)) (* r (sin angle) 0+i))) angles))) ; Convert the polar coordinates to complex numbers

;; This function calculates the sigma value for the Aberth method
(define (calculate-sigma approximations k)
  (let ((zk (list-ref approximations k))) ; Get the k-th approximation
    (foldl (lambda (zj sigma)             ; Fold over the approximations
             (if (not (= zj zk))          ; If zj is not equal to zk
                 (+ sigma (/ 1 (- zk zj))) ; Add 1/(zk - zj) to sigma
                 sigma))                  ; Otherwise, keep sigma unchanged
           0.0                            ; Start with sigma = 0
           (remove (lambda (x) (equal? x (list-ref approximations k))) approximations)))) ; Remove zk from the list of approximations

;; This function implements the Aberth method for finding polynomial roots
(define (aberth-method coefficients (epsilon 0.0001)) ; Default epsilon is 0.0001
  (let* ((approximations (get-approximations coefficients)) ; Get initial approximations
         (derivative-coefficients (polynomial-derivative-coefficients coefficients)) ; Compute the derivative coefficients
         (max-iterations 100))            ; Set the maximum number of iterations

    (define (iterate approximations iteration) ; Define a recursive function to iterate the Aberth method
      (if (>= iteration max-iterations)   ; If the maximum number of iterations is reached
          approximations                  ; Return the current approximations
          (let ((offsets
                 (for/list ((zk approximations) (k (in-naturals))) ; Loop over the approximations
                   (let* ((frac (frac-val coefficients derivative-coefficients zk)) ; Compute f(zk)/f'(zk)
                          (sigma (calculate-sigma approximations k))) ; Compute sigma_k
                     (/ frac (- 1 (* frac sigma))))))) ; Compute the correction term
            (let ((new-approximations (map - approximations offsets))) ; Update the approximations
              (define (check-convergence approximations) ; Define a function to check if the roots have converged
                (foldl (lambda (val roots-converge) ; Fold over the approximations
                         (if (not (is-close-to-zero? (poly-val (reverse coefficients) val) epsilon)) ; Check if the polynomial value is close to zero
                             #f                 ; If not, return false
                             roots-converge))   ; Otherwise, keep checking
                       #t                       ; Start with roots-converge = true
                       approximations))         ; Fold over the list of approximations
              (if (check-convergence new-approximations) ; If the roots have converged
                  new-approximations           ; Return the new approximations
                  (iterate new-approximations (+ iteration 1))))))) ; Otherwise, continue iterating

    (iterate approximations 0)))                ; Start the iteration with the initial approximations and iteration 0

;; Read polynomial coefficients from a file
(define coefficients (read-coefficients "poly_coeff(997).txt"))
;; Measure the start time
(define start (current-inexact-milliseconds))
;; Find the roots using the Aberth method
(define result (aberth-method coefficients))
;; Measure the end time
(define end (current-inexact-milliseconds))
;; Print the time taken for the computation
;;(printf "Time: ~a ms\n" (- end start))
;; Uncomment the following lines to print the roots and their rounded versions
;; (printf "Roots: ~a\n" result)
(printf "Rounded Roots: ~a\n" (round-complex result))
;; (printf "len ~a\n" (length result))
