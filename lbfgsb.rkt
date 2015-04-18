#lang racket/base

(require ffi/unsafe
	 ffi/vector
	 ffi/cvector
	 racket/list
	 racket/match
	 math/flonum)

(define lbfgsb (ffi-lib "lbfgsb.so"))

(define setulb
  (get-ffi-obj 'setulb_ lbfgsb
	       (_fun (n : (_ptr i _int))
		     (m : (_ptr i _int))
		     (x : _pointer)
		     (l : (_f64vector i))
		     (u : (_f64vector i))
		     (nbd : (_s32vector i))
		     (f : (_ptr i _double))
		     (g : _pointer)
		     (factr : (_ptr i _double))
		     (pgtol : (_ptr i _double))
		     (wa : (_f64vector i))
		     (iwa : (_s32vector i))
		     (task : _bytes)
		     (iprint : (_ptr i _int))
		     (csave : _bytes)
		     (lsave : (_s32vector i))
		     (isave : (_s32vector i))
		     (dsave : (_f64vector i))
		     (task_len : _int = (bytes-length task))
		     (csave_len : _int = (bytes-length csave))
		     -> _void)))

(define (create-wa n m)
  (list->f64vector (make-list (+ (* 2 m n)
				 (* 11 m m)
				 (* 5 n)
				 (* 8 n)) 0)))

(define (create-iwa n m)
  (list->s32vector (make-list (* 3 n) 0)))

(define (approximate-fprime f epsilon)
  (lambda (x [f0 #f])
    (define finit (if (number? f0)
		      f0
		      (f x)))
    (define n (flvector-length x))
    (define d (make-flvector n 0.0))
    (define xtemp (flvector-copy x))
    
    (for/flvector ([i (in-naturals)]
		   [v x])
      (flvector-set! xtemp i (fl+ v epsilon))
      (begin0
       (/ (- (f xtemp) finit) epsilon)
       (flvector-set! xtemp i v)))))

(define (minimize fun x0
		  #:jac [jac #f]
		  #:bounds [bounds #f]
		  #:ftol [ftol 2.22e-9]
		  #:gtol [gtol 1e-5]
		  #:eps [epsilon 1e-8]
		  #:iprint [iprint -1]
		  #:maxcor [maxcor 10]
		  #:maxfun [maxfun 15000]
		  #:maxiter [maxiter 15000])
  (define x (for/flvector ([x x0]) x))
  (define n (flvector-length x))
  (define m maxcor)
  (define l (list->f64vector (make-list n 0)))
  (define u (list->f64vector (make-list n 0)))
  (define nbd (list->s32vector (make-list n 0)))

  (cond
   [(equal? bounds #f) (void)]
   [(not (= (length bounds) n)) (error "bounds length mismatch")]
   [(= (length bounds) n)
    (for ([bound bounds]
	  [i (in-naturals)])
      (match bound
	[(list (? number? lo) (? number? hi))
	 (f64vector-set! l i lo)
	 (f64vector-set! u i hi)
	 (s32vector-set! nbd i 2)]
	[(list _ (? number? hi))
	 (f64vector-set! u i hi)
	 (s32vector-set! nbd i 3)]
	[(list (? number? lo) _)
	 (f64vector-set! l i lo)
	 (s32vector-set! nbd i 1)]))])

  (printf "~a ~a ~a ~a ~a ~a~n"
	  (flvector->list x)
	  n
	  m
	  (f64vector->list l)
	  (f64vector->list u)
	  (s32vector->list nbd))
  
  (define factr (/ ftol epsilon.0))
  (define pgtol gtol)
  
  (define wa (create-wa n m))
  (define iwa (create-iwa n m))
  (define task (make-bytes 60 32))
  (bytes-copy! task 0 #"START")
  (define csave (make-bytes 60 32))
  (define lsave (make-s32vector 4))
  (define isave (make-s32vector 44))
  (define dsave (make-f64vector 29))

  (define gfun
    (if (procedure? jac)
	jac
	(approximate-fprime fun epsilon)))
  
  (define f 0.0)
  (define g (list->flvector (make-list n 0)))

  (define-values
    (funcalls iters)
    (let/ec break
      (for/fold
       ([funcalls 0]
	[iter 0])
       ([ignore (in-naturals)])
       #:break (or (> funcalls maxfun)
		   (> iter maxiter))
       
       (setulb n m
	       (flvector->cpointer x) l u nbd
	       f (flvector->cpointer g)
	       factr pgtol
	       wa iwa
	       task iprint
	       csave lsave isave dsave)
       
       (cond
	[(equal? #"FG" (subbytes task 0 2))
	 (begin
	   (set! f (fun x))
	   (set! g (gfun x f))
	   (values (add1 funcalls)
		   (add1 iter)))]
	[(equal? #"NEW_X" (subbytes task 0 5))
	 (values funcalls
		 (add1 iter))]
	[else (break funcalls iter)]))))
  
  (cond
   [(equal? #"CONV" (subbytes task 0 4))
    x]
   [(equal? #"ABNO" (subbytes task 0 4))
    (error "abnormal termination")]
   [(equal? #"ERROR" (subbytes task 0 5))
    (error "error in input parameters")]
   [else
    (error (bytes->string/utf-8 task))]))

(flvector->list
 (minimize (lambda (x) (flsin (flvector-ref x 0)))
	   (vector 0.5)
	   #:bounds '((-0.5 #f)))
 )
