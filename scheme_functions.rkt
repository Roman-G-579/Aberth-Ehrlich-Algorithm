#lang racket

(define (polynomial-eval coeffs x)
  (let loop ((coeffs coeffs) (result 0))
    (if (null? coeffs)
        result
        (loop (cdr coeffs) (+ (* result x) (car coeffs))))))

(define (polynomial-derivative coeffs)
  (let loop ((coeffs (cdr coeffs)) (result '()) (n (sub1 (length coeffs))))
    (if (null? coeffs)
        (reverse result)
        (loop (cdr coeffs) (cons (* (car coeffs) n) result) (sub1 n)))))

(define (random-approximations coeffs)
  (define n (sub1 (length coeffs)))
  (define radius (expt (/ (abs (car (reverse coeffs))) (abs (car coeffs))) (/ 1 n)))
  (define angles (map (lambda (k) (* k (/ (* 2 pi) n))) (range n)))
  (map (lambda (angle) (make-rectangular (* radius (cos angle)) (* radius (sin angle)))) angles))

(define (parse-list str)
  (map string->number (string-split str ",")))

(define (format-list lst)
  (string-join (map number->string lst) ","))

(define (format-complex-list lst)
  (string-join (map (lambda (c) (string-append (number->string (real-part c)) "," (number->string (imag-part c)))) lst) ","))

(define (parse-complex str)
  (let* ((parts (string-split str ",")))
    (make-rectangular (string->number (first parts)) (string->number (second parts)))))

(define (parse-complex-list str)
  (map parse-complex (string-split str ",")))

(define (main func args output-file)
  (define result
    (case (string->symbol func)
      ((polynomial-eval) (let ((coeffs (parse-list (first args)))
                               (x (string->number (second args))))
                           (polynomial-eval coeffs x)))
      ((polynomial-derivative) (let ((coeffs (parse-list (first args))))
                                 (polynomial-derivative coeffs)))
      ((random-approximations) (let ((coeffs (parse-list (first args))))
                                 (random-approximations coeffs)))
      (else (error "Unknown function"))))
  (call-with-output-file output-file
    (lambda (out)
      (if (list? result)
          (if (complex? (first result))
              (fprintf out "~a" (format-complex-list result))
              (fprintf out "~a" (format-list result)))
          (fprintf out "~a" result)))
    #:exists 'replace))

;; Get command-line arguments and call main
(let ([args (cdr (current-command-line-arguments))])
  (main (list-ref args 0) (list (list-ref args 1)) (list-ref args 2)))
