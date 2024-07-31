#lang racket

(define (round-complex complex-nums)
  (map (lambda (num)
         (let ((real (real-part num))
               (imag (imag-part num)))
           (+ (/ (round (* real 1000)) 1000)
              (* (/ (round (* imag 1000)) 1000) 0+i))))
       complex-nums))

(define (read-coefficients filename)
  (with-input-from-file filename
    (lambda ()
      (let loop ((lines '()))
        (let ((line (read-line)))
          (if (eof-object? line)
              (reverse lines)
              (let ((num (string->number (string-trim line))))
                (if (number? num)
                    (loop (cons num lines))
                    (begin
                      (printf "Error: Invalid number format in line: ~a\n" line)
                      (loop lines))))))))))

(define (is-close-to-zero? value epsilon)
  (and (< (abs (real-part value)) epsilon)
       (< (abs (imag-part value)) epsilon)))
(define pi 3.141592653589793)

 (define poly-val (λ(p x)
         (if (null? p)
             0
             (+ (car p) (* x (poly-val (cdr p) x)))
         )
 ))

(define (frac-val p q x)
  (if (<= (magnitude x) 1)
      (/ (poly-val (reverse p) x) (poly-val (reverse q) x))
      (* x (/ (poly-val p (/ 1 x))
              (poly-val q (/ 1 x))))))

(define (polynomial-derivative-coefficients coefficients)
  (let* ((n (length coefficients))
         (derivative-coeffs '()))
    (if (= n 1)
        '(0)
        (for/list ((i (in-range 0 (- n 1))))
          (* (- n 1 i) (list-ref coefficients i))))))

(define (get-approximations coefficients)
  (let* ((deg (- (length coefficients) 1))
         (pi-val pi)
         (p-0 (list-ref coefficients deg))
         (p-n (car coefficients))
         (r (magnitude (expt (/ p-0 p-n) (/ 1 deg))))
         (angles (for/list ((x (in-range 0 deg))) (* x (/ (* 2 pi-val) deg)))))
    (map (lambda (angle) (+ (* r (cos angle)) (* r (sin angle) 0+i))) angles)))

(define (calculate-sigma approximations k)
  (let ((zk (list-ref approximations k)))
    (foldl (lambda (zj sigma)
             (if (not (= zj zk))
                 (+ sigma (/ 1 (- zk zj)))
                 sigma))
           0.0
           (remove (lambda (x) (equal? x (list-ref approximations k))) approximations))))

(define (aberth-method coefficients (epsilon 0.0001))
  (let* ((approximations (get-approximations coefficients))
         (derivative-coefficients (polynomial-derivative-coefficients coefficients))
         (max-iterations 100))

    (define (iterate approximations iteration)
      (if (>= iteration max-iterations)
          approximations
          (let ((offsets
                 (for/list ((zk approximations) (k (in-naturals)))
                   (let* ((frac (frac-val coefficients derivative-coefficients zk))
                          (sigma (calculate-sigma approximations k)))
                     (/ frac (- 1 (* frac sigma)))))))
            (let ((new-approximations (map - approximations offsets)))
              (define (check-convergence approximations)
                (foldl (lambda (val roots-converge)
                         (if (not (is-close-to-zero? (poly-val (reverse coefficients) val) epsilon))
                             #f
                             roots-converge))
                       #t
                       approximations))
              (if (check-convergence new-approximations)
                  new-approximations
                  (iterate new-approximations (+ iteration 1)))))))

    (iterate approximations 0)))


; Example usage
(define coefficients (read-coefficients "test2.txt"))
(define start (current-inexact-milliseconds))
(define result (aberth-method coefficients))
(define end (current-inexact-milliseconds))
(printf "Time: ~a ms\n" (- end start))
(printf "Roots: ~a\n" result)
(printf "Rounded Roots: ~a\n" (round-complex result))
(printf "len ~a\n" (length result))