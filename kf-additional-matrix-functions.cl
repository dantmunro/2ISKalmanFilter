;; -=Begin Copyright Notice=-
;; copyright (c) 2002-2015 2Is Inc, Walpole, MA  All rights reserved.
;;
;; All rights reserved.
;;
;; Restricted Rights Legend
;; ------------------------
;; Use, duplication, and disclosure of the software, data and information
;; contained herein by any agency, department or entity of the U.S.
;; Government are subject to restrictions of Restricted Rights for
;; Commercial Software developed at private expense as specified in FAR
;; 52.227-19 or DOD FAR Supplement 252.227-7013 (c) (1) (ii), as
;; applicable.
;; -=End Copyright Notice=-

;; A series of in-house matrix functions comprising the 2IS Statistics library.
;; Not authored by Daniel Munro, unlike all of the preceding code

(in-package :kf)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;            MACROS

(defmacro 2dbl (number)
  `(if (typep ,number 'ratio) ,number (coerce ,number 'double-float)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; initialization functions ;;;

;;; functions create new matrix
(defun make-matrix-w-normal-random-vals (m n)
  (let ((matrix (make-array (list m n))))
    (init-array-w-normal-random-vals matrix)
    matrix))
	
;;; functions modify array in place
(defun init-array-w-normal-random-vals (a)
  (init-array-w-fun a #'(lambda () (normal-random-variable))))

(defun init-array-w-uniform-random-vals (a)
  (init-array-w-fun a #'(lambda ()
			  (uniform-random-variable
			   :uniform-input (random 1d0)))))

(defun init-array-w-fun (a f)
  (loop for i from 0 to (1- (array-total-size a)) do
	(setf (row-major-aref a i) (funcall f)))
  a)

;;; unary matrix functions ;;;

;; multiply a matrix by a number (scalar)
(defun matrix-times-scalar (matrix scalar)
  (apply-fun-to-array #'(lambda (n) (* scalar n)) matrix))

;; apply unary function f to each element of array a creating a new array new-a.
(defun apply-fun-to-array (f a)
  (let ((new-a (make-array (array-dimensions a))))
    (loop for i from 0 to (1- (array-total-size a))
	do (setf (row-major-aref new-a i) (funcall f (row-major-aref a i))))
    new-a))

;; make a copy of a matrix
(defun copy-matrix (m)
  (copy-array m))

;; make a copy of an array
(defun copy-array (a)
  (let ((new-a (make-array (array-dimensions a) :element-type (array-element-type a))))
    (loop for i from 0 to (1- (array-total-size a))
	do (setf (row-major-aref new-a i) (row-major-aref a i)))
    new-a))
	
;;; binary matrix functions ;;;

;; add 2 matrices together
;; this function assumes that the dimensions of the two matrices are the same
(defun matrix-add (matrix1 matrix2)
  (apply-fun-to-array-pair #'+ matrix1 matrix2))

;; subtract 2 matrices
;; this function assumes that the dimensions of the two matrices are the same
(defun matrix-subtraction (matrix1 matrix2)
  (apply-fun-to-array-pair #'- matrix1 matrix2))

;; apply binary function f to each element of array a creating a new array new-a.
;; Both arrays are assumed to have the same dimensions.
(defun apply-fun-to-array-pair (f a1 a2)
  (let ((new-a (make-array (array-dimensions a1))))
    (loop for i from 0 to (1- (array-total-size a1)) do
	(setf (row-major-aref new-a i)
	  (funcall f (row-major-aref a1 i) (row-major-aref a2 i))))
    new-a))
	
;;; reduce matrix functions ;;;

; multiply 2 matrices; # of num-cols of first matrix must equal # of row of second matrix
(defun matrix-multiply (matrix1 matrix2)
  (let* ((m (array-dimension matrix1 0))
         (n (array-dimension matrix1 1))
         (p (array-dimension matrix2 1))
         (c (make-array (list m p))))
    (loop for i from 0 to (1- m)
        do
          (loop for k from 0 to (1- p)
              do
                (setf (aref c i k) (loop for j from 0 to (1- n)
                                       summing
                                         (* (aref matrix1 i j) (aref matrix2 j k))))))
    c))

; transpose a matrix
(defun matrix-transpose (matrix)
  (let* ((m (array-dimension matrix 0))
         (n (array-dimension matrix 1))
         (c (make-array (list n m))))
    (loop for i from 0 to (1- m)
        do
          (loop for j from 0 to (1- n)
              do
                (setf (aref c j i) (aref matrix i j))))
    c))

; determinant of a 2 x 2 matrix
(defun 2x2matrix-det (matrix)
  (if (and (= (array-dimension matrix 0) 2)
           (= (array-dimension matrix 1) 2))
      (- (* (2dbl (aref matrix 0 0)) (2dbl(aref matrix 1 1))) 
         (* (2dbl (aref matrix 0 1)) (2dbl (aref matrix 1 0))))))

; determinant of a 3 x 3 matrix
(defun 3x3matrix-det (matrix)
  (if (and (= (array-dimension matrix 0) (array-dimension matrix 1)) 
           (= (array-dimension matrix 0) 3))
      (+ (* (aref matrix 0 0) (aref matrix 1 1) (aref matrix 2 2)) (* (aref matrix 0 2) (aref matrix 1 0) (aref matrix 2 1)) (* (aref matrix 0 1) (aref matrix 1 2) (aref matrix 2 0))
         (- (* (aref matrix 0 0) (aref matrix 1 2) (aref matrix 2 1))) (- (* (aref matrix 0 1) (aref matrix 1 0) (aref matrix 2 2)))
         (- (* (aref matrix 0 2) (aref matrix 1 1) (aref matrix 2 0))))))

; must be used for a square matrix
(defun minor-matrix (matrix rownumber colnumber)
  (if (> (array-dimension matrix 0) colnumber)
      (make-array (list (1- (array-dimension matrix 0)) 
                        (1- (array-dimension matrix 1))) 
                  :initial-contents (loop for i from 0 to (1- (array-dimension matrix 0))     
                                        when (not (= i rownumber))
                                        collect
                                          (loop for j from 0 to (1- (array-dimension matrix 1))
                                              when (not (= j colnumber))
                                              collect
                                                (aref matrix i j))))))

; calculates the determinant of a matrix; do not use this function for a non-square matrix
(defun determinant (matrix)
  (if* (= 1 (array-dimension matrix 0))
     then
          (aref matrix 0 0)
   elseif (= 2 (array-dimension matrix 0))
     then
          (2x2matrix-det matrix)
     else
          (loop for j from 0 to (1- (array-dimension matrix 1))
              when (not (= (array-dimension matrix 0) 2))
              sum
                (* (expt -1 j) (2dbl (aref matrix 0 j)) (determinant (minor-matrix matrix 0 j))))))
				
(defun matrix-inverse-1x1 (matrix)
  (let* ((value (aref matrix 0 0))
         (new-value (/ 1 value)))
    (make-array '(1 1) :initial-contents (list (list new-value)))))
				
(defun matrix-inverse-2x2 (matrix)
  (let ((det (determinant matrix))
        (inv (make-array '(2 2))))
    (setf (aref inv 0 0) (aref matrix 1 1)
      (aref inv 1 1) (aref matrix 0 0)
      (aref inv 0 1) (- (aref matrix 0 1))
      (aref inv 1 0) (- (aref matrix 1 0)))
    (matrix-times-scalar inv (/ det))))

; can only calculate the inverse of a square matrix
; uses Gauss-Jordan elimination to determine the inverse
(defun matrix-inverse (matrix)
  (let* ((num-rows (array-dimension matrix 0))
         (num-cols (array-dimension matrix 1)))
    (cond ((/= num-rows num-cols)
           (error "cannot invert a non-square matrix."))
          ((= num-rows 1)
           (matrix-inverse-1x1 matrix))
          ((= num-rows 2)
           (matrix-inverse-2x2 matrix))
          (t
           (let ((xnum-cols (+ num-rows num-cols))
                 (vector (make-array num-rows))
                 (inverse (make-array (list num-rows num-cols))))
             (loop for r from 0 to (1- num-rows)
                 as v = (make-array xnum-cols :initial-element 0)
                 do (setf (aref v (+ r num-rows)) 1)
                   (loop for c from 0 to (1- num-cols)
                       do
                         (setf (aref v c) (aref matrix r c))
                         (setf (aref vector r) v)))
             (loop for r from 0 to (1- num-rows)
                 do
                   (loop for i from r to (1- num-rows)
                       with imax = -1
                       with vmax = -1
                       as value = (abs (aref (aref vector i) r))
                       when (> value vmax)
                       do (setq vmax value)
                         (setq imax i)
                       finally 
                         (if (zerop vmax) (return-from matrix-inverse ()))
                         (if* (/= r imax)
                            then
                                 (let ((temp-change (aref vector r)))
                                   (setf (aref vector r) (aref vector imax))
                                   (setf (aref vector imax) temp-change))))
                   (let ((pivot-row (aref vector r)))
                     (loop for i from (1+ r) to (1- xnum-cols)
                         with divisor = (aref pivot-row r)
                         do 
                           (setf (aref pivot-row i) (/ (aref pivot-row i) divisor)))
                     (loop for i from (1+ r) to (1- num-rows)
                         as row = (aref vector i)
                         do 
                           (loop for j from (1+ r) to (1- xnum-cols)
                               with divisor = (aref row r)
                               do 
                                 (setf (aref row j) (- (aref row j) (* divisor (aref pivot-row j))))))))
             (loop for r from (1- num-rows) downto 1
                 as pivot-row = (aref vector r)
                 do
                   (loop for i from (1- r) downto 0
                       as row = (aref vector i)
                       do
                         (loop for j from num-rows to (1- xnum-cols)
                             with divisor = (aref row r)
                             do 
                               (setf (aref row j) (- (aref row j) (* divisor (aref pivot-row j)))))))
             (loop for r from 0 to (1- num-rows)
                 as v = (aref vector r)
                 do
                   (loop for c from 0 to (1- num-cols)
                       do 
                         (setf (aref inverse r c) (aref v (+ c num-cols)))))
             inverse)))))
			 
(defun uniform-random-variable (&key (uniform-input (random 1.0))
                                     (min 0) (max 1.0) (delta (- max min)))
  (+ min (* uniform-input delta)))
 
(defun normal-random-variable (&optional (mu 0) (sigma 1))
  ;box-muller transform
  (let ((U1 (random 1.0))
        (U2 (random 1.0)))
    (if (zerop U1) (normal-random-variable) ;a remote possibility
      (+ mu
         (* sigma
            (sqrt (* -2 (log U1)))
            (cos (* 2 pi U2))))))) ;cosine is allegedly expensive.
			 
(defun identity-matrix (n)
  (let ((x (make-array (list n n) :initial-element 0)))
    (loop for i from 0 to (- n 1)
        do (setf (aref x i i) 1))
    x))
	
(defun print-matrix (matrix &optional name)
  (let ((n (array-dimension matrix 0))
        (m (array-dimension matrix 1)))
    (if name (format t "~&~A:~%" name))
    (loop for i from 0 to (- n 1)
        do (loop for j from 0 to (- m 1)
               do (format t "~vt~a" (* j 5) (aref matrix i j))
               finally (terpri)))))