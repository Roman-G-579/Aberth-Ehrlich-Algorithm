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

(define (cos x)
  (define (cos-helper x term n cos-x)
    (if (> n 25)
        cos-x
        (cos-helper x
                    (* term (- (/ (* x x) (* (- (* 2 n) 1) (* 2 n)))))
                    (+ n 1)
                    (+ cos-x term))))
  (cos-helper x 1 1 1))

(define (sin x)
  (define (sin-helper x term n sin-x)
    (if (> n 25)
        sin-x
        (sin-helper x
                    (* term (- (/ (* x x) (* (* 2 n) (+ (* 2 n) 1)))))
                    (+ n 1)
                    (+ sin-x term))))
  (sin-helper x x 1 x))

(define pi 3.141592653589793)

(define (poly-val p x)
  (let loop ((p p) (x x) (s 0))
    (if (null? p)
        s
        (loop (cdr p) x (+ (car p) (* x s))))))

(define (frac-val p q x)
  (if (<= (magnitude x) 1)
      (/ (poly-val p x) (poly-val q x))
      (* x (/ (poly-val (reverse p) (/ 1 x))
              (poly-val (reverse q) (/ 1 x))))))

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

(define (aberth-method coefficients (epsilon 0.0001))
  (let* ((approximations (get-approximations coefficients))
         (derivative-coefficients (polynomial-derivative-coefficients coefficients)))
    (define (iterate approximations)
      (let ((offsets
             (for/list ((zk approximations) (k (in-naturals)))
               (let* ((frac (frac-val coefficients derivative-coefficients zk))
                      (sigma (for/sum ((zj approximations) (j (in-naturals)) #:when (not (= k j)))
                               (/ 1 (- zk zj)))))
                 (/ frac (- 1 (* frac sigma)))))))
        (let ((new-approximations
               (map - approximations offsets)))
          (if (< (apply max (map magnitude offsets)) epsilon)
              new-approximations
              (iterate new-approximations)))))
    (iterate approximations)))

; Example usage
(define coefficients (read-coefficients "poly_coeff(997).txt"))
(define start (current-inexact-milliseconds))
(define result (aberth-method coefficients))
(define end (current-inexact-milliseconds))
(printf "Time: ~a ms\n" (- end start))
;(printf "Roots: ~a\n" result)
(printf "Rounded Roots: ~a\n" (round-complex result))
