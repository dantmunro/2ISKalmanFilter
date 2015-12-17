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

(in-package :kf)

; Only intended for 2d-arrays and vectors of 2d arrays. Directs output to output-to-csv-mat
; or output-to-csv-vec
(defun output-to-csv (item item-name
						&key (out-file "/debug-lisp.csv"))
	(if (vectorp item) (output-to-csv-vec item item-name :out-file out-file) 
						(output-to-csv-mat item item-name :out-file out-file))
)

; Outputs a 2D array to a CSV with a header (note that it is appended, not supeseded)
(defun output-to-csv-mat (mat mat-name 
							&key (out-file "/debug-lisp.csv"))
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
							&key (out-file "/debug-lisp.csv"))
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
