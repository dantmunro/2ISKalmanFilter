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

(in-package :kf)

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