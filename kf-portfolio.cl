;; Single runnable bash script containing all modules necessary for 2IS Kalman 
;; Filter.  This script is intended to serve as an entirely self-contained 
;; Kalman Filter/Latent Linear Dynamical System (LDS) implementation.

;; This project was designed for Allegro Common Lisp implementation but 
;; the script is for use with CLISP

;#!/usr/local/bin/clisp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
;;
;; Build $Revision: 15 $
;;
;; -=End Copyright Notice=-

;; A series of previously unimplemented statistical and matrix functions, along with a series of 
;; shorthand functions for already existing ones in the statistics package

;; Constants
(defconstant ROWS 0) ; Indices in an array-dimensions matrix represnting rows, columns respectively
(defconstant COLUMNS 1)
(defconstant 1D 1)
(defconstant 2D 2)
(defconstant KF-PARAMS '(transition-matrix emission-matrix transition-mean prior-mean
							emission-mean transition-covariance prior-covariance emission-covariance))
(defconstant DEFAULT-TOL 1.e-5)
(defconstant DEFAULT-MAXITS 100)

;; Some Matrix Algebra Shorthand	
(defun %+% (matrix1 matrix2 &rest other-matrices)
	(if (null other-matrices)
		(matrix-add matrix1 matrix2)
		(reduce #'matrix-add (append (list matrix1 matrix2) other-matrices))))

(defun %-% (matrix1 &optional matrix2)
	(if (null matrix2) 
		(%-% (zeros matrix1) matrix1)
		(matrix-subtraction matrix1 matrix2)))

; N.B. Using this operator on >2 operands has occasionally returned slightly different
; results than performing it on only 2.  As many statements as possible have been stripped out
; which involve the use of multiple %*% operators where only one would theoretically do.
; However, keep this in mind in case any unforseen resulting errors arise.
(defun %*% (matrix1 matrix2 &rest other-matrices) 
	(if (null other-matrices)
		(matrix-multiply matrix1 matrix2)
		(reduce #'matrix-multiply (append (list matrix1 matrix2) other-matrices))))
		
(defun %.*% (var1 var2) ; Multiply constant by matrix (in either order)
	(cond ((and (arrayp var1) (numberp var2)) (matrix-times-scalar var1 var2))
		((and (arrayp var2) (numberp var1)) (matrix-times-scalar var2 var1))))
	
(defun inv (matrix)
	(matrix-inverse matrix))
	
(defun %/% (matrix1 matrix2)
	(%*% matrix1 (inv matrix2)))
	
(defun %//% (matrix1 matrix2)
	(%*% (inv matrix1) matrix2))
	
(defun tr (matrix)
	(matrix-transpose matrix))
	
(defun printmat (matrix)
	(terpri) (print-matrix matrix))
	
(defun eye (n)
	(identity-matrix n))
	
(defun det (matrix)
	(determinant matrix))
	
(defun copy (matrix)
	(copy-matrix matrix))
	
(defun randmat (r c)
	(make-matrix-w-normal-random-vals r c))

;; Additional statistical/matrix operations
	
; Make a matrix of zeros with the following bases:
; 1. If dim1 and dim2 are both integers, make a matrix with dim1 rows and dim2 cols
; 2. If dim1, dim2 and dim3 are all integers, make a dim1xdim2xdim3 matrix
; 3. If only dim1 is non-null, and it is an integer, make a dim1xdim1
; 4. If dim1 is an array, make an array of zeros with its dimensions
; 5. If dim1 is a list, make an array of zeros with dimensions matching the lists'
; 	first two elements
(defun zeros (dim1 &optional dim2 dim3)
  (cond 
	((and (integerp dim1) (integerp dim2) (null dim3)) 
		(make-array (list dim1 dim2) :initial-element 0d0))
	((and (integerp dim1) (integerp dim2) (integerp dim3))
		(make-array (list dim1 dim2 dim3) :initial-element 0d0))
	((integerp dim1)
		(make-array (list dim1 dim1) :initial-element 0d0))
	((arrayp dim1)
		(make-array (array-dimensions dim1) :initial-element 0d0))
	((listp dim1)
		(make-array (subseq dim1 0 2) :initial-element 0d0))))

; Analagous function for ones
(defun ones (dim1 &optional dim2 dim3)
  (cond 
	((and (integerp dim1) (integerp dim2) (null dim3)) 
		(make-array (list dim1 dim2) :initial-element 1d0))
	((and (integerp dim1) (integerp dim2) (integerp dim3))
		(make-array (list dim1 dim2 dim3) :initial-element 1d0))
	((integerp dim1)
		(make-array (list dim1 dim1) :initial-element 1d0))
	((arrayp dim1)
		(make-array (array-dimensions dim1) :initial-element 1d0))
	((listp dim1)
		(make-array (subseq dim1 0 2) :initial-element 1d0))))
					  
; 2D array subindexing (nondestructive). Takes an array and up to two 2-tuple lists denoting the 
; horizontal and vertical ranges of the matrix/array to subindex. (if either the row-range or 
; col-range is left blank, the selection will include all rows or columns, respectively).
; Returns the respective subarray in 2D-array form.
(defun 2d-aref-range (arr &optional row-range col-range)
	(cond
		((null row-range) (setf row-range (list 0 (1- (array-dimension arr ROWS)))))
		((integerp row-range) (setf row-range (list row-range row-range))))
	(cond
		((null col-range) (setf col-range (list 0 (1- (array-dimension arr COLUMNS)))))
		((integerp col-range) (setf col-range (list col-range col-range))))
	(let ((2d-list (loop for i from (first row-range) to (second row-range) 
		collect (loop for j from (first col-range) to (second col-range)
				collect (aref arr i j)))))
		(make-array (list (length 2d-list)
                    (length (first 2d-list)))
					:initial-contents 2d-list)))
	
; Analagously modifies a subarray of a 2d array (destructive)	
(defun 2d-array-modify (arr to-replace-with &optional row-range col-range)
	(cond
		((null row-range) (setf row-range (list 0 (1- (array-dimension arr ROWS)))))
		((integerp row-range) (setf row-range (list row-range row-range))))
	(cond
		((null col-range) (setf col-range (list 0 (1- (array-dimension arr COLUMNS)))))
		((integerp col-range) (setf col-range (list col-range col-range))))
	(loop for i from (first row-range) to (second row-range)
		for k below (array-dimension to-replace-with ROWS)
		do (loop for j from (first col-range) to (second col-range)
			for l below (array-dimension to-replace-with COLUMNS)
			do (setf (aref arr i j) (aref to-replace-with k l))))
	arr)

; Converts a vector of column matrices to a single matrix
(defun concat-columns-to-mat (vec)
	(let ((concat (zeros (array-total-size (aref vec 0)) (length vec))))
		(loop for j below (array-dimension concat COLUMNS)
			do
			(loop for i below (array-dimension concat ROWS)
				do
				(setf (aref concat i j) (aref (aref vec j) i 0))
			)
		)
		concat
	)
)
	
(defun is-vector (arr) ;Checks if an array is a vector, or a 2D array where one dimension is 1
	(or (= (length (array-dimensions arr)) 1) 
		(loop for dim in (array-dimensions arr) thereis (= dim 1))))

(defun sqrt-1d-array (vec) ;Takes a row or column 2D array and takes elementwise square root
	(let* ((tosqrt vec))
		(cond ((= 1 (array-dimension tosqrt ROWS)) 
				(loop for i below (array-dimension tosqrt COLUMNS) 
					do (setf (aref tosqrt 0 i) (sqrt (aref tosqrt 0 i))) finally (return tosqrt)))
				((= 1 (array-dimension tosqrt COLUMNS))
				(loop for i below (array-dimension tosqrt ROWS) 
					do (setf (aref tosqrt i 0) (sqrt (aref tosqrt i 0))) finally (return tosqrt)))
				(t (error "Vector must be 1-D")))
	)
)
	
(defun norm (x)
  (let ((len (car (array-dimensions x))))
    (sqrt (loop for i from 0 to (1- len) sum (expt (aref x i 0) 2)))))

; Givens Rotation - Characteristic operation of the Givens QR Decomposition Algorithm
(defun givens-rotation (a b)
	(let ((c) (s) (r))
		(if (zerop b)
			(progn
				(setf c 1)
				(setf s 0)
			)
			(progn
				(if (> (abs b) (abs a))
					(progn
						(setf r (/ a b))
						(setf s (/ 1 (sqrt (+ 1 (expt r 2)))))
						(setf c (* s r))
					)
					(progn
						(setf r (/ b a))
						(setf c (/ 1 (sqrt (+ 1 (expt r 2)))))
						(setf s (* c r))
					)
				)
			)
		)
		(values c s)
	)
)

; QR Decomposition with Givens Rotations 
(defun qr-decomp-givens (A)
	(let* ((m (array-dimension A ROWS)) (n (array-dimension A COLUMNS))
		(Q (eye m)) (R (copy A)))
		(loop for j below n
			do 
			(loop for i from (1- m) downto (1+ j)
				for G = (eye m)
				do
				(multiple-value-bind (c s) (givens-rotation (aref R (1- i) j) (aref R i j))
					(setf (aref G (1- i) (1- i)) c)
					(setf (aref G i (1- i)) s)
					(setf (aref G (1- i) i) (- s))
					(setf (aref G i i) c)
					(setf R (%*% (tr G) R))
					(setf Q (%*% Q G))
				)
			)
		)
		(loop for i below m
			do
			(loop for j below n
				if (< j i)
				do (setf (aref R i j) 0d0)
			)
		)
		(values Q R)
	)
)

; Get upper triangular matrix of a matrix
(defun triu (A k)
	(let ((triu-A (zeros `,(array-dimensions A))))
		(loop for i below (array-dimension A ROWS)
			do
			(loop for j below (array-dimension A COLUMNS)
				if (>= j (+ k i)) 
				do (setf (aref triu-A i j) (aref A i j))
			)
		)
		triu-A
	)
)

; Concatenate entire matrix into a column
(defun column-concat (A)
	(let ((concat-vector (make-array (array-total-size A) :fill-pointer 0))
			(concat-matrix (zeros (array-total-size A) 1)))
		(loop for j below (array-dimension A COLUMNS)
			do
			(loop for i below (array-dimension A ROWS)
				do
				(vector-push (aref A i j) concat-vector)
			)
		)
		(loop for i below (length concat-vector)
			do
			(setf (aref concat-matrix i 0) (aref concat-vector i))
		)
		concat-matrix
	)
)

; Determines the dimensionality of arr, and does one of the following:
; 1. If arr is effectively 1D, creates a diagonal matrix out of arr
; 2. If arr is 2D, returns a column of arr's diagonal elements
(defun diag (arr)
	(cond
		((or (= (length (array-dimensions arr)) 1D)
				(loop for dim in (array-dimensions arr) thereis (= dim 1))) (vec-to-diag arr))
		((= (length (array-dimensions arr)) 2D) (diag-to-vec arr))
		(t (error "Error: Argument must be array"))
	)
)

; Performs the first of diag's two operations
(defun vec-to-diag (arr)
	(let ((diag-mat (zeros (array-total-size arr) (array-total-size arr))))
		(loop for i below (array-total-size arr)
			do
			(setf (aref diag-mat i i) 
				(if (= (length (array-dimensions arr)) 2D) ; Array is 2D
					(if 
						(= (first (array-dimensions arr)) 1) (aref arr 0 i) ; Array is 2D row
						(aref arr i 0) ; Array is 2D column
					)
					(aref arr i) ; Array is 1D
				)
			)
		)
		diag-mat
	)
)

; Performs the second of diag's two operations
(defun diag-to-vec (matrix)
	(let* ((diag-length  (apply #'min (array-dimensions matrix)))
			(diag-vec (zeros diag-length 1)))
		(loop for i below diag-length
			do
			(setf (aref diag-vec i 0) (aref matrix i i))
		)
		diag-vec
	)
)

; Performs Paul Godfrey's Singular Value Decomposition, which decomposes a matrix a into
; a matrix of a's eigenvectors (u), a diagonal matrix of a's eigenvalues (s), and the application
; of the PCA (Principal Component Analysis) (v)
(defun svdsim (a &key (tol (* double-float-epsilon 1024)) (debug-p nil))
	(let* ((size-a (array-dimensions a)) 
		(loop-max (* 100 (reduce #'max size-a)))
		(loop-count 0) (u (eye (first size-a))) (s (tr a))
		(v (eye (second size-a))) (err most-positive-long-float)
		(q (zeros `,(array-dimensions a))))
		(loop while (and (> err tol) (< loop-count loop-max))
			do
			(multiple-value-setq (q s) (qr-decomp-givens (tr s)))
			when debug-p do 
				(loop for mat in (list q s)
					for mat-name in (list "q" "s")
					do (printmat-with-label mat mat-name)
				)
			do (setf u (%*% u q))
			when debug-p do (printmat-with-label u "u")
			do (multiple-value-setq (q s) (qr-decomp-givens (tr s)))
			when debug-p do 
				(loop for mat in (list q s)
					for mat-name in (list "q" "s")
					do (printmat-with-label mat mat-name)
				)
			(setf v (%*% v q))
			when debug-p do (printmat-with-label v "v")
			do
			(let* ((e-mat (triu s 1)) 
				(E (norm (column-concat e-mat)))
				(F (norm (tr (diag s)))))
				(when debug-p (progn (printmat-with-label e-mat "e")
									(print "E=") (terpri) (print E) (terpri)
									(print "F=") (terpri) (print F) (terpri)))
				(when (zerop F) (setf F 1))
				(setf err (/ E F))
				(when debug-p (print "err=") (terpri) (print err) (terpri))
				(incf loop-count)
			)
		)
		(when debug-p 
			(loop for mat in (list u s v)
				for mat-name in (list "u" "s" "v")
				do (printmat-with-label mat mat-name)
			)
		)
		(let ((ss (diag s)) (s (zeros (reduce #'max `,size-a)))) 
			(loop for n below (array-dimension ss ROWS)
				for ssn = (aref ss n 0)
				do
				(setf (aref s n n) (abs ssn))
				if (< 0 ssn)
					do
					(loop for i below (array-dimension u ROWS)
						do
						(setf (aref u i n) (- (aref u i n)))
					)
			)
			(values (%-% u) s v)
		)
	)
)

(defun mean (xi) (/ (reduce #'+ xi) (length xi)))

(defun variance (xi) (/ (reduce #'+ (mapcar (lambda (x) (expt (- x (mean xi)) 2)) xi)) (1- (length xi))))

(defun covariance (xi xj) 
	(/ (reduce #'+ (mapcar (lambda (x y) (* (- x (mean xi)) (- y (mean xj)))) xi xj)) (1- (length xi))))
	
(defun covariance-matrix (matrix)
	(let ((S (zeros (array-dimension matrix COLUMNS) (array-dimension matrix COLUMNS))))
		(loop for i below (array-dimension matrix COLUMNS)
			for xi = 
				(loop for l below (array-total-size 
					(2d-aref-range matrix nil i)) collect (aref (2d-aref-range matrix nil i) l 0))
			do
			(setf (aref S i i) (variance xi))
			(loop for j from (1+ i) to (1- (array-dimension matrix COLUMNS))
				for xj = 
				(loop for l below (array-total-size 
					(2d-aref-range matrix nil j)) collect (aref (2d-aref-range matrix nil j) l 0)) 
				do
				(setf (aref S i j) (covariance xi xj))
				(setf (aref S j i) (aref S i j))
			)
		)
		S
	)
)

; Draws a vector from a multivariate normal distribution
(defun mvrandn (mu sigma)
	(let ((D (array-dimension mu ROWS)))
		(multiple-value-bind (U S V) (svdsim sigma)
			(let ((r (randmat D 1)))
				(setf r (%+% (%*% U (%*% (diag (sqrt-1d-array (diag S))) r)) mu))
				r
			)
		)
	)
)

; Nondestructively concatenates two matrices horizontally
(defun append-matrices-horiz (mat1 mat2)
	(let ((r1 (array-dimension mat1 ROWS)) (c1 (array-dimension mat1 COLUMNS))
			(r2 (array-dimension mat2 ROWS)) (c2 (array-dimension mat2 COLUMNS)))
		(assert (= r1 r2))
		(let ((newmat (make-array (list r1 (+ c1 c2)))))
			(loop for i below r1
				do (loop for j below c1
					do (setf (aref newmat i j) (aref mat1 i j))
				)
				(loop for j below c2
					do (setf (aref newmat i (+ c1 j)) (aref mat2 i j))
				)
			)
			newmat
		)
	)
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
;;
;; Build $Revision: 2 $
;;
;; -=End Copyright Notice=-

;; A handful of methods created to output 2d-arrays and vectors of 2d arrays to CSV files, 
;; for debugging purposes (perhaps to compare with analagous output from MATLAB/Octave)

; Only intended for 2d-arrays and vectors of 2d arrays. Directs output to output-to-csv-mat
; or output-to-csv-vec
(defun output-to-csv (item item-name
						&key (out-file "~/plottable_data/debug-lisp.csv"))
	(if (vectorp item) (output-to-csv-vec item item-name :out-file out-file) 
						(output-to-csv-mat item item-name :out-file out-file))
)

; Outputs a 2D array to a CSV with a header (note that it is appended, not supeseded)
(defun output-to-csv-mat (mat mat-name 
							&key (out-file "~/plottable_data/debug-lisp.csv"))
	(with-open-file (str out-file
					 :direction :output
					 :if-exists :append
					 :if-does-not-exist :create)
		(format str "~a~%" mat-name)
		(output-mat str mat)
	)
)

; Outputs a vector of 2D arrays to a CSV with a header (note that it is appended, not supeseded)
(defun output-to-csv-vec (vec vec-name
							&key (out-file "~/plottable_data/debug-lisp.csv"))
	(with-open-file (str out-file
					 :direction :output
					 :if-exists :append
					 :if-does-not-exist :create)
		(format str "~a~%" vec-name)
		(loop for a across vec do
			(if (numberp a) ;Assumes every vector contains only a number or 2d array
				(progn (format str "~,,,,,,'eE," a) (terpri str))
				(output-mat str a)
			)
		)
	)
)

; Outputs just a matrix, with no header
(defun output-mat (str mat)
	(loop for i below (array-dimension mat ROWS)
		do
		(loop for j below (array-dimension mat COLUMNS)
			do (format str "~,,,,,,'eE," (aref mat i j))
		)
		(terpri str)
	)
)

; Delete a file, only if it already exists
(defun delete-if-exists (path) (when (not (null (probe-file path))) (delete-file path)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;;; initialization functions ;;;

;;; functions create new matrix
(defun make-matrix-w-normal-random-vals (m n)
  (let ((matrix (make-array (list m n))))
    (init-array-w-normal-random-vals matrix)
    matrix))

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

; multiply 2 matrices; # of cols of first matrix must equal # of row of second matrix
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
  (let* ((rows (array-dimension matrix 0))
         (cols (array-dimension matrix 1)))
    (cond ((/= rows cols)
           (error "cannot invert a non-square matrix."))
          ((= rows 1)
           (matrix-inverse-1x1 matrix))
          ((= rows 2)
           (matrix-inverse-2x2 matrix))
          (t
           (let ((xcols (+ rows cols))
                 (vector (make-array rows))
                 (inverse (make-array (list rows cols))))
             (loop for r from 0 to (1- rows)
                 as v = (make-array xcols :initial-element 0)
                 do (setf (aref v (+ r rows)) 1)
                   (loop for c from 0 to (1- cols)
                       do
                         (setf (aref v c) (aref matrix r c))
                         (setf (aref vector r) v)))
             (loop for r from 0 to (1- rows)
                 do
                   (loop for i from r to (1- rows)
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
                     (loop for i from (1+ r) to (1- xcols)
                         with divisor = (aref pivot-row r)
                         do 
                           (setf (aref pivot-row i) (/ (aref pivot-row i) divisor)))
                     (loop for i from (1+ r) to (1- rows)
                         as row = (aref vector i)
                         do 
                           (loop for j from (1+ r) to (1- xcols)
                               with divisor = (aref row r)
                               do 
                                 (setf (aref row j) (- (aref row j) (* divisor (aref pivot-row j))))))))
             (loop for r from (1- rows) downto 1
                 as pivot-row = (aref vector r)
                 do
                   (loop for i from (1- r) downto 0
                       as row = (aref vector i)
                       do
                         (loop for j from rows to (1- xcols)
                             with divisor = (aref row r)
                             do 
                               (setf (aref row j) (- (aref row j) (* divisor (aref pivot-row j)))))))
             (loop for r from 0 to (1- rows)
                 as v = (aref vector r)
                 do
                   (loop for c from 0 to (1- cols)
                       do 
                         (setf (aref inverse r c) (aref v (+ c cols)))))
             inverse)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;; Additional Matrix/statistics functions.  Also not authored by Daniel Munro.

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
			   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
			   
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
;;
;; Build $Revision: 37 $
;;
;; -=End Copyright Notice=-



;;; Kalman Filter (Latent Linear Dynamical System/LDS) Forward/Backward and Training Algorithms
;;; Includes some demos and several auxiliary matrix/statistical functions
;;; Developers: Daniel Munro

(in-package :kf)

(defun printmat-with-label (mat mat-name) (format t "~a=" mat-name) (terpri) (printmat mat) (terpri))

; Kalman Filter constructor.  This is the standard order for field referencing, meaning that for 
; readability, it is requested that all future developers - when given the choice - reference
; these fields in the following order when performing similar operations on a large number of them
(defclass kalman-filter ()
	((observations				:accessor v :initarg :v :initform nil)
	(dim-filter 				:accessor H :initarg :H :initform 0)
	(dim-observations			:accessor V :initarg :V :initform 0)
	(num-observations			:accessor T :initarg :T :initform 0)
	(transition-matrix			:accessor A :initarg :A :initform nil)
	(emission-matrix			:accessor B :initarg :B :initform nil)
	(transition-mean			:accessor meanH :initarg :meanH :initform nil)
	(prior-mean					:accessor meanP :initarg :meanP :initform nil)
	(emission-mean				:accessor meanV :initarg :meanV :initform nil)
	(transition-covariance		:accessor CovH :initarg :CovH :initform nil)
	(prior-covariance			:accessor CovP :initarg :CovP :initform nil)
	(emission-covariance		:accessor CovV :initarg :CovV :initform nil)
	(filtered-means				:accessor f :initarg :f :initform nil)
	(filtered-covariances		:accessor F :initarg :F :initform nil)
	(smoothed-means				:accessor g :initarg :g :initform nil)
	(smoothed-covariances		:accessor G :initarg :G :initform nil)
	(smoothed-2-step-posteriors	:accessor Gp :initarg :Gp :initform nil)
	(log-predicted-prob			:accessor logpgvgv :initarg :logpgvgv :initform nil)
	(sequence-log-likelihood	:accessor L :initarg :L :initform 0))
)

; Kalman Filter constructor, part 2.  Note that all output arrays (f,F,g,G,Gp) are *fixed-size*.
; In other words, the Kalman Filter can only operate on a single static observation matrix
; Since 3D arrays in Lisp have a counterintuitive orientation (each member 2D array is considered
; its own row rather than a tier) I have decided to use vectors of 2D arrays for f, F and logpgvgv,
; rather than 3D arrays 
(defmethod initialize-instance :after ((filter kalman-filter) &key)
	(with-accessors ((v v) (V V) (T T) (H H) (f f) (F F) (g g) (G G) (Gp Gp) (logpgvgv logpgvgv)) filter
		(setf V (array-dimension v ROWS))
		(setf T (array-dimension v COLUMNS))
		(assert (and (> H 0) (> V 0) (> T 0)))
		
		; The following five lines cannot be factored due to accessor restrictions
		(setf f (make-array T :initial-element (zeros H 1)))
		(setf F (make-array T :initial-element (zeros H H)))
		(setf g (make-array T :initial-element (zeros H 1)))
		(setf G (make-array T :initial-element (zeros H H)))
		(setf Gp (make-array (1- T) :initial-element (zeros H H)))
		(setf logpgvgv (make-array T :initial-element 0))
	)
)

; Used for the learning algorithm.  Can take a Kalman Filter (only meant to hold params
; A,B,meanH,meanP,meanV,CovH,CovP,CovV), a tolerance and/or max iterations, all three of which can 
; replace the defaults which would be supplised into a running of the learning algorithm
(defstruct (learning-params (:conc-name params-)) filter tolerance maxits)

; Forward Algorithm - calculates the core Predict/Update Kalman Filter operations and returns
; filtered state estimates (means and covariances)
; Slots used - v,A,B,CovH,CovV,meanH,meanV,CovP,meanP
; Slots returned - modifications of f (p(h(i|v(1:i))), F (p(h(i)|v(1:i))),
; L (log likelihood of sequence log p(v(1:T)))
(defmethod forward ((filter kalman-filter) &key (debug-p nil))
	(with-accessors ((T T) (f f) (F F) (logpgvgv logpgvgv) (L L)) filter
		(when debug-p (print "Starting forward"))
		(loop for i below T do (forward-update filter i :debug-p debug-p))
		(when debug-p (print "Finished forward"))
		(setf L (reduce #'+ (logpgvgv filter)))
		(when debug-p 
			(progn 
				(loop for mat in (list f F)
					for mat-name in (list "f" "F")
					do (format t "~a=~%" mat-name) (mapcar #'printmat (coerce mat 'list)) (terpri)
				)
				(print "L = ") (print L) (terpri)
			)
		) 
	)
)

; A single unit of the forward algorithm - performs it for a given state i
; Slots used - f,F,v,A,B,CovH,CovV,meanH,meanV
; (In this case, v is a single observation, so it is represented by its index (obs-ind) in the actual
; v matrix)
; Slots returned - modifications of f, F and "logpgvgv" ( = log p(v(i+1)|v(1:i)))
(defmethod forward-update ((filter kalman-filter) obs-ind &key (debug-p nil))
	(with-accessors ((A A) (CovH CovH) (CovV CovV)) filter
		(when (is-vector A) (setf A (diag A)))
		(when (is-vector CovH) (setf CovH (diag CovH)))
		(when (is-vector CovV) (setf CovV (diag CovV)))
	)
	(let* ( ;I avoided a with-accessors statement here to avoid confusion and let statement conflicts
		(f (if (zerop obs-ind) (zeros (H filter) 1) (aref (f filter) (1- obs-ind))))
		(F (if (zerop obs-ind) (zeros (H filter) (H filter)) (aref (F filter) (1- obs-ind))))
		(v (2d-aref-range (v filter) nil obs-ind))
		(meanH (if (zerop obs-ind) (meanP filter) (meanH filter)))
		(CovH (if (zerop obs-ind) (CovP filter) (CovH filter)))
		(muh (%+% (%*% (A filter) f) meanH))
		(muv (%+% (%*% (B filter) muh) (meanV filter)))
		(Shh (%+% (%*% (A filter) F (tr (A filter))) CovH))
		(Svv (%+% (%*% (B filter) Shh (tr (B filter))) (CovV filter)))
		(Svh (%*% (B filter) Shh))
		(del (%-% v muv))
		(invSvvdel (%//% Svv del))
		(fnew (%+% muh (%*% (tr Svh) invSvvdel)))
		(K (%/% (%*% Shh (tr (B filter))) Svv))
		(tmp (%-% (eye (array-dimension (A filter) ROWS)) (%*% K (B filter))))
		(Fnew (%+% (%*% tmp Shh (tr tmp)) (%*% K (CovV filter) (tr K))))) ; Joseph's Form
		(setf (aref (f filter) obs-ind) fnew)
		(setf (aref (F filter) obs-ind) Fnew)
		(setf (aref (logpgvgv filter) obs-ind) (+ (aref (%.*% -.5 (%*% (tr del) invSvvdel)) 0 0)
													(* -.5 (log (det Svv)))
													(* (array-dimension v ROWS) -.5 (log (* 2 pi)))))
		(when debug-p (print "logpgvgv =") (loop for l across (logpgvgv filter) do (print l)) (terpri))
	)
)

; Backward Algorithm - performs the Rauch-Tung-Striebel Smoother and returns "smoothed" state estimates
; that are more refined than the filtered estimates.
; Slots used - v,A,B,f,F,CovH,meanH
; Slots returned - modifications of g, G
(defmethod backward ((filter kalman-filter) &key (debug-p nil))
	(with-accessors ((T T) (f f) (F F) (g g) (G G) (Gp Gp)) filter
		(when debug-p (print "Starting backward"))
		(setf (aref g (1- T)) (aref f (1- T)))
		(setf (aref G (1- T)) (aref F (1- T)))
		(loop for i from (- T 2) downto 0 do (backward-update filter i))
		(when debug-p
			(print "Finished backward")
			(loop for mat in (list g G Gp)
				for mat-name in (list "g" "G" "Gp")
				do (format t "~a=~%" mat-name) (mapcar #'printmat (coerce mat 'list)) (terpri)
			)
		)
	)
)

; A single unit of the backward algorithm - performs it for a given state i
; Slots used - g,G,f,F,A,CovH,meanH
; Slots returned - modifications of g, G, Gp (smoothed mean, cov, cross moment <h_t h_{i+1}|v(1:T)>)
(defmethod backward-update ((filter kalman-filter) obs-ind)
	(with-accessors ((A A) (CovH CovH)) filter
		(when (is-vector A) (setf A (diag A)))
		(when (is-vector CovH) (setf CovH (diag CovH)))
	)
	(let* ( ;I avoided a with-accessors statement here to avoid confusion and let statement conflicts
		(g (aref (g filter) (1+ obs-ind)))
		(G (aref (G filter) (1+ obs-ind)))
		(f (aref (f filter) obs-ind))
		(F (aref (F filter) obs-ind))
		(muh (%+% (%*% (A filter) f) (meanH filter)))
		(Shtptp (%+% (%*% (A filter) F (tr (A filter))) (CovH filter)))
		(Shtpt (%*% (A filter) F))
		(leftA (%/% (tr Shtpt) Shtptp))
		(leftS (%-% F (%*% leftA Shtpt)))
		(leftm (%-% f (%*% leftA muh)))
		(gnew (%+% (%*% leftA g) leftm))
		(Gnew (%+% (%*% leftA G (tr leftA)) leftS))
		(Gnew (%.*% .5 (%+% Gnew (tr Gnew))))
		(Gpnew (%+% (%*% leftA G) (%*% gnew (tr g)))))
		(setf (aref (g filter) obs-ind) gnew)
		(setf (aref (G filter) obs-ind) Gnew)
		(setf (aref (Gp filter) obs-ind) Gpnew)
	)
)

; Just an encapsulation of the full forward-backward process
(defmethod smooth ((filter kalman-filter) &key (debug-p nil) (persist-p nil))
	(with-accessors ((T T) (v v) (f f) (F F) (g g) (G G) (Gp Gp)) filter
		;Refresh T in case synthetic observations added
		(setf T (array-dimension v COLUMNS)) 
		(forward filter :debug-p debug-p)
		(backward filter :debug-p debug-p)
		(when persist-p 
			(loop for mat in (list f F g G Gp)
				for mat-name in (list "f" "F" "g" "G" "Gp")
				do (output-to-csv-vec mat mat-name)
			)
		)
	)
)

; Learn the parameters A,B,meanH,meanP,meanV,CovH,CovP,CovV through the EM algorithm
; (where the E-step is equivalent to smoothing/forward-backward)
; N.B. opt-filter = "varargin" or "varargin{1}" in Matlab/Octave code
(defun learn (v H &key (opt-filter nil) (debug-p nil) (persist-p nil))
	;Delete any files that would likely be outputted to if "persist-p" was turned on 
	;(this will have to be changed when the output directory does - it is recommended that
	;anyone who codes extensions to the Kalman filter keeps such a directory)
	(when persist-p
		(loop for i in '("debug-lisp" "original-lisp-data" "reconstructed-lisp-data") do
			(delete-if-exists (make-pathname 
									:directory '(:absolute "Users" "dmunro" "Downloads" "plottable_data") 
									:name (format nil "~a" i) :type "csv"))
		) 
	)
	(let* ((init (make-instance 'kalman-filter :v v :H H)) 
			(V (array-dimension v ROWS)) (T (array-dimension v COLUMNS)))
		(with-accessors ((A A) (B B) (meanH meanH) (meanP meanP) (meanV meanV) 
				(CovH CovH) (CovP CovP) (CovV CovV)) init
			(setf A (eye H)) ;Giving default starting values to each of the parameters
			(multiple-value-bind (usvd ssvd vsvd) (svdsim v) ;Get the principal components thru SVD
				(when debug-p 
					(loop for mat in (list usvd ssvd vsvd)
						for mat-name in (list "usvd" "ssvd" "vsvd")
						do (printmat-with-label mat mat-name)
					)
				)
				(setf B (2d-aref-range usvd nil (list 0 (1- H))))
				(setf CovH (eye H))
				(setf CovP (eye H))
				(let ((vv (%*% (2d-aref-range usvd nil (list 0 (1- H))) 
						(%*% (2d-aref-range ssvd (list 0 (1- H)) (list 0 (1- H))) 
							(2d-aref-range vsvd (list 0 (1- H)) nil)))))
					(when debug-p
						(loop for mat in (list vv (%-% v vv))
							for mat-name in (list "vv" "(%-% v vv)")
							do (printmat-with-label mat mat-name)
						)
					)
					(setf CovV (covariance-matrix (tr (%-% v vv)))) 
					(setf meanP (zeros H 1))
					(setf meanH (zeros H 1))
					(setf meanV (zeros V 1))
				)
			)
		)
		;Create a "learning-params" object to hold any non-default parameters (for the filter's
		;parameters, the max iterations of EM, or the convergence tolerance).  Though it would have
		;been more expedient code-wise to give "tolerance" and "maxits" defaults in the defstruct
		;statement, the need to choose between the defaults for such values or manually inputted ones
		;makes this unattainable.
		(let* ((opts (make-learning-params  
				:filter ;Merge any hardcoded filter parameters from opt-filter with the 
						;default ones from init
					(loop for slot in KF-PARAMS
						with new-filter = init
						when (and (not (null opt-filter))
								(not (null (slot-value (params-filter opt-filter) slot))))
						do (setf (slot-value init slot) (slot-value (params-filter opt-filter) slot))
						finally (return new-filter)
					)
				:tolerance (if (or (null opt-filter) (null (params-tolerance opt-filter))) DEFAULT-TOL 
								(params-tolerance opt-filter))
				:maxits (if (or (null opt-filter) (null (params-maxits opt-filter))) DEFAULT-MAXITS 
							(params-maxits opt-filter))))  
				(loglikold (make-array (params-maxits opts) :fill-pointer 0)) 
				(diagCovH) (diagCovV) (diagCovP) (diagA) 
				(loglik (make-array (params-maxits opts) :fill-pointer 0)))
			(vector-push most-negative-long-float loglikold)
			(with-accessors ((A A) (B B) (meanH meanH) (meanP meanP) (meanV meanV) 
					(CovH CovH) (CovP CovP) (CovV CovV)) init 
				(loop for mat in (list A B CovH CovP CovV meanH meanP meanV)
					for mat-name in (list "A" "B" "CovH" "CovP" "CovV" "meanH" "meanP" "meanV")
					when debug-p do (printmat-with-label mat mat-name) 
					when persist-p do (output-to-csv-mat mat mat-name)
				)
			)
			(with-accessors ((A A) (B B) (meanH meanH) (meanV meanV) (meanP meanP) (CovV CovV)
					(CovH CovH) (CovP CovP) (f f) (F F) (g g) (G G) (Gp Gp) (logpgvgv logpgvgv) (L L))
					(params-filter opts)
				(when (is-vector CovH) (setf diagCovH t))
				(when (is-vector CovV) (setf diagCovV t))
				(when (is-vector CovP) (setf diagCovP t))
				(when (is-vector A) (setf diagA t))
				(loop for loop1 below (params-maxits opts)
					do
					(when debug-p (progn (print "ITERATION #") (print loop1) (terpri)))
					(smooth (params-filter opts) :debug-p nil :persist-p nil) ;smoothing - E-step
					(vector-push L loglik)
					when debug-p 
						do (print "loglik=") (loop for i across loglik do (print i)) (terpri) (terpri)
					when (and (> loop1 0) (< (aref loglik loop1) (aref loglik (1- loop1))))
						do (print "Warning: Log Likelihood Decreased")
					when (loop for i across loglikold 
							always (< (- (aref loglik loop1) i) (params-tolerance opts)))
						do (return)
					unless (loop for i below (array-total-size meanP) always (zerop (aref meanP i 0))) 
						do (setf meanP (aref g 0))
					do ;M-Step PHASE 1: Recalculate CovP and meanP (prior distribution)
					(setf CovP (%+% 
								(aref G 0) (%*% (aref g 0) (tr (aref g 0)))
								(%-% (%*% (aref g 0) (tr meanP)))
								(%-% (%*% meanP (tr (aref g 0))))		
								(%*% meanP (tr meanP)))
					)
					(setf CovP (%.*% .5 (%+% CovP (tr CovP))))
					when diagCovP do (setf CovP (diag (diag CovP)))
					do
					(loop for mat in (list meanP CovP) 
						for mat-name in (list "meanP" "CovP")
						when debug-p do (printmat-with-label mat mat-name)
						when persist-p do (output-to-csv-mat mat mat-name)
					)
					;M-Step PHASE 2: Recalculate A and B (transition/emission matrices)
					(let ((HH (zeros H H)) (HHp (zeros H H)) (VH (zeros V H)))
						(loop for t1 below (- T 1)
							do
							(setf HH (%+% HH (aref G t1) (%*% (aref g t1) (tr (aref g t1)))))
							(setf HHp (%+% HHp (tr (aref Gp t1)) (%-% (%*%  meanH (tr (aref g t1))))))
							(setf VH (%+% VH 
										(%*% (%-% (2d-aref-range v nil t1) meanV) (tr (aref g t1))))
							)
							finally 
							(loop for mat in (list HH HHp VH)
								for mat-name in (list "HH" "HHp" "VH")
								when debug-p do (printmat-with-label mat mat-name)
								when persist-p do (output-to-csv-mat mat mat-name)
							)
						)
						(setf A (%/% HHp HH))
						(when diagA (setf A (diag (diag A))))
						(setf HH (%+% HH (aref G (1- T)) (%*% (aref g (1- T)) (tr (aref g (1- T))))))
						(setf VH (%+% VH (%*% 
											(%-% (2d-aref-range v nil (1- T)) meanV)
											(tr (aref g (1- T)))))
						)
						(setf B (%/% VH HH))
						(loop for mat in (list A B)
							for mat-name in (list "A" "B")
							when debug-p do (printmat-with-label mat mat-name)
							when persist-p do (output-to-csv-mat mat mat-name)
						)
					)
					;M-Step PHASE 3: Recalculate CovV and meanV (emission/observation distribution)
					(let ((CovVtmp (zeros V V)) (meanVtmp (zeros V 1)))
						(loop for t1 below T
							for dv = (%-% (2d-aref-range v nil t1) meanV)
							do
							(setf CovVtmp (%+% 
											CovVtmp 
											(%*% dv (tr dv)) 
											(%-% (%*% dv (%*% (tr (aref g t1)) (tr B))))
											(%-% (%*% B (%*% (aref g t1) (tr dv))))
											(%*% B 
												(%*% 
													(%+% (aref G t1) (%*% (aref g t1) (tr (aref g t1))))
													(tr B))))
							)
							(setf meanVtmp 
								(%+% meanVtmp (2d-aref-range v nil t1) (%-% (%*% B (aref g t1)))))
						)
						(when (not (loop for i below (array-total-size meanV) 
							always (zerop (aref meanV i 0))))  
							(setf meanV (%.*% (/ 1 T) meanVtmp))
						)
						(setf CovV (%.*% (/ 1 T) CovVtmp))
						(setf CovV (%.*% (%+% CovV (tr CovV)) .5))
						(when diagCovV (setf CovV (diag (diag CovV))))
					)
					(loop for mat in (list meanV CovV)
						for mat-name in (list "meanV" "CovV")
						when debug-p do (printmat-with-label mat mat-name)
						when persist-p do (output-to-csv-mat mat mat-name)
					)
					;M-Step PHASE 4: Recalculate CovH and meanH (transition distribution)
					(let ((CovHtmp (zeros H H)) (meanHtmp (zeros H 1)))
						(loop for t1 below (- T 1)
							do
							(setf CovHtmp 
								(%+% 
									CovHtmp (aref G (1+ t1)) (%*% (aref g (1+ t1)) (tr (aref g (1+ t1))))
									(%-% (%*% A (aref Gp t1))) (%-% (%*% (tr (aref Gp t1)) (tr A)))
									(%*% A 
										(%*% 
											(%+% (aref G t1) (%*% (aref g t1) (tr (aref g t1)))) 
											(tr A)))
								)
							)
							(setf CovHtmp
								(%+%
									CovHtmp (%-% (%*% (aref g (1+ t1)) (tr meanH))) 
									(%-% (%*% meanH (tr (aref g (1+ t1))))) (%*% meanH (tr meanH))
									(%*% meanH (%*% (tr (aref g t1)) (tr A)))
									(%*% A (aref g t1) (tr meanH))
								)
							)
							(setf meanHtmp (%+% meanHtmp (aref g (1+ t1)) (%-% (%*% A (aref g t1)))))
						)
						(setf CovH (%.*% CovHtmp (/ 1 (1- T))))
						(setf CovH (%.*% .5 (%+% CovH (tr CovH))))
						(when diagCovH (setf CovH (diag CovH)))
						(when (not (loop for i below (array-total-size meanH) 
							always (zerop (aref meanH i 0)))) (setf meanH (%.*% (/ 1 (1- T)) meanHtmp)))
					)
					(loop for mat in (list meanH CovH)
						for mat-name in (list "meanH" "CovH")
						when debug-p do (printmat-with-label mat mat-name)
						when persist-p do (output-to-csv-mat mat mat-name)
					)
					(setf loglikold loglik)
					(when persist-p (output-to-csv-vec loglik "loglik"))
				)
			)
		)
		init
	)
)


; Use learned parameters and a partial sequence of observations to incrementally extrapolate future 
; hidden state estimates using synthetic observations - artificially generated through the emission 
; equation (v(i) = B*g(i) + N(meanV,CovV))
; N.B. predict is a nondestructive method
(defmethod predict ((filter kalman-filter) num-new-preds)
	(loop for i below num-new-preds 
		with v = (v filter)
		for newprior = (%*% (A filter) (aref (g filter) (1- (length (g filter)))))
		for newobs = (%+% (%*% (B filter) newprior) 
							(mvrandn (zeros (array-dimension v ROWS) 1) (CovV filter)))
		do
		; Since the size of the state estimate vectors (f,F,g,G,Gp) are fixed, append the new
		; observation to the observation matrix, and build a new Kalman Filter around it, passing
		; in the same learned parameters
		(setf v (append-matrices-horiz v newobs))
		(setf filter (make-instance 'kalman-filter :H (H filter) :v v :A (A filter) 
			:B (B filter) :CovH (CovH filter) :CovP (CovP filter) :CovV (CovV filter) 
			:meanH (meanH filter) :meanP (meanP filter) :meanV (meanV filter)))
		(smooth filter)
		finally (return filter)
	)
)

; Equivalent to demoLDSlearn.m
(defun demo-LDS-learn ()
	(let* ((T 50)
			(lat (make-array `(1 ,T) 
				:initial-contents (list (loop for i from 1 to T collect (sin (* .2 i))))))
			(B (make-matrix-w-normal-random-vals 2 1))
			(v (%*% B lat))
			(v (%+% v (%.*% .1 (make-matrix-w-normal-random-vals 
						(array-dimension v ROWS) (array-dimension v COLUMNS)))))
			(H 2)
			(opts (make-learning-params
					:maxits 100
					:filter (make-instance 'kalman-filter :H H :v v :CovV (ones 2 1))))
		)
		(let* ((filter (learn v H :opt-filter opts :debug-p nil :persist-p nil)))
			(smooth filter)
			(let ((vv (%*% (B filter) (concat-columns-to-mat (g filter)))))
				(with-open-file (str 
								"~/plottable_data/original-lisp-data.csv"
								 :direction :output
								 :if-exists :supersede
								 :if-does-not-exist :create)
					(loop for i below (array-dimension v ROWS)
						do
						(loop for j below (array-dimension v COLUMNS)
							do (format str "~,,,,,,'eE," (aref v i j))
						)
						(terpri str)
					)
				)
				(with-open-file (str 
								"~/plottable_data/reconstructed-lisp-data.csv"
								 :direction :output
								 :if-exists :supersede
								 :if-does-not-exist :create)
					(loop for i below (array-dimension vv ROWS)
						do
						(loop for j below (array-dimension vv COLUMNS)
							do (format str "~,,,,,,'eE," (aref vv i j))
						)
						(terpri str)
					)
				)
				(loop for mat in (list v vv)
					for mat-name in (list "v" "vv")
					do (printmat-with-label mat mat-name)
				)
			)
		)
	)
)

; Equivalent to demoLDStracking.m
(defun demo-LDS-tracking ()
	(let* ((T 400) (Delta .1) (A (eye 6)) (B #2A((0d0 1d0 0d0 0d0 0d0 0d0) (0d0 0d0 0d0 1d0 0d0 0d0))))
		(loop for ind in (list '(0 4) '(1 0) '(2 5) '(3 2))
			do
			(setf (aref A (car ind) (cadr ind)) Delta)
		)
		(let ((h (zeros 6 T)) (sigV 50d0) (sigH .00001) (v (zeros 2 T))) ;FIXME Convert sigH to double
			(setf (aref h 1 0) (random 1d0))
			(setf (aref h 3 0) (random 1d0))
			(setf (aref h 0 0) (* 15 (random 1d0)))
			(setf (aref h 2 0) (* 15 (random 1d0)))
			(setf (aref h 4 0) (random 1d0))
			(setf (aref h 5 0) (- (random 1d0)))
			(2d-array-modify v (%+% (%*% B (2d-aref-range h nil 0)) (%.*% sigV (randmat 2 1))) nil 0)
			(loop for t1 from 1 to (1- T)
				do
				(2d-array-modify h 
					(%+% (%*% A (2d-aref-range h nil (1- t1))) (%.*% sigH (randmat 6 1))) nil t1)
				(2d-array-modify v
					(%+% (%*% B (2d-aref-range h nil t1)) (%.*% sigV (randmat 2 1))) nil t1)
			)
			(let* ((CovH (%.*% (* sigH sigH) (eye 6))) (CovV (%.*% (* sigV sigV) (eye 2)))
					(CovP (eye 6)) (meanP (zeros 6 1)) (meanH (zeros 6 1)) (meanV (zeros 2 1))
					(filter (make-instance 'kalman-filter :H (array-dimension h ROWS) :v v
								:A A :B B :meanH meanH :meanP meanP :meanV meanV :CovH CovH 
								:CovP CovP :CovV CovV)))
				(smooth filter)
				(with-accessors ((A A) (B B) (meanH meanH) (meanV meanV) (meanP meanP) (CovV CovV)
					(CovH CovH) (CovP CovP) (f f) (F F) (g g) (G G) (Gp Gp) (logpgvgv logpgvgv) (L L))
						filter
					(loop for mat in (list f F g G Gp v h)
							for mat-name in (list "dum1" "dum2" "mean_post" "cov_post" "dum3" "v" "h")
						do (output-to-csv mat mat-name
							:out-file 
								(format nil "~/plottable_data/~a.csv" mat-name))
					)
				)
			)
		)
	)
)

(defun demo-LDS-prediction (obs-num num-new-preds)
	(let* ((T obs-num)
			(lat (make-array `(1 ,T) 
				:initial-contents (list (loop for i from 1 to T collect (sin (* .2 i))))))
			(B (make-matrix-w-normal-random-vals 2 1))
			(v (%*% B lat))
			(v (%+% v (%.*% .1 (make-matrix-w-normal-random-vals 
						(array-dimension v ROWS) (array-dimension v COLUMNS)))))
			(H 2)
			(opts (make-learning-params
					:maxits 100
					:filter (make-instance 'kalman-filter :H H :v v :CovV (ones 2 1))))
		)
		(let* ((filter (learn v H :opt-filter opts :debug-p nil :persist-p t)))
			(smooth filter)
			(setf filter (predict filter num-new-preds))
			(with-accessors ((v v) (A A) (B B) (meanH meanH) (meanV meanV) (meanP meanP) (CovV CovV)
					(CovH CovH) (CovP CovP) (f f) (F F) (g g) (G G) (Gp Gp) (logpgvgv logpgvgv) (L L))
					  filter
				(let ((vv (%*% B (concat-columns-to-mat g))))
					(loop for mat in (list v vv)
						for filename in (list "original-lisp-data1" "reconstructed-lisp-data1")
						do
						(with-open-file (str (format nil 
												"~/plottable_data/~a.csv"
												filename)
										 :direction :output
										 :if-exists :supersede
										 :if-does-not-exist :create)
							(loop for i below (array-dimension mat ROWS)
								do
								(loop for j below (array-dimension mat COLUMNS)
									do (format str "~,,,,,,'eE," (aref mat i j))
								)
								(terpri str)
							)
						)
					)
					(loop for mat in (list v vv)
						for mat-name in (list "v" "vv")
						do (printmat-with-label mat mat-name)
					)
				)
			)
		)
	)
)

;;;UNCOMMENT TO RUN DEMOS
(demo-LDS-prediction 30 5)
;(demo-LDS-learn)
;(demo-LDS-tracking)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;