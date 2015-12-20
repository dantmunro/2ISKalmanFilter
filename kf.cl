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

;; All packages intended for use with Allegro Common Lisp

;;; Kalman Filter (Latent Linear Dynamical System/LDS) Forward/Backward and Training Algorithms
;;; Includes some demos and several auxiliary matrix/statistical functions
;;; Developers: Daniel Munro

(in-package :kf)
#|
;; Constants
(defconstant ROWS 0) ; Indices in an array-dimensions matrix represnting rows, columns respectively
(defconstant COLUMNS 1)
(defconstant 1D 1)
(defconstant 2D 2)
(defconstant KF-PARAMS '(transition-matrix emission-matrix transition-mean prior-mean
							emission-mean transition-covariance prior-covariance emission-covariance))
(defconstant DEFAULT-TOL 1.e-5)
(defconstant DEFAULT-MAXITS 100)
|#

(defun printmat-with-label (mat mat-name) (format t "~a=" mat-name) (terpri) (printmat mat) (terpri))

; Kalman Filter constructor.  This is the standard order for field referencing, meaning that for 
; readability, it is requested that all future developers - when given the choice - reference
; these fields in the following order when performing similar operations on a large number of them
(defclass kalman-filter ()
	((observations				:accessor v :initarg :v :initform nil)
	(dim-filter 				:accessor H :initarg :H :initform 0)
	(dim-observations			:accessor V :initarg :V :initform 0)
	(num-observations			:accessor T1 :initarg :T1 :initform 0)
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
	(with-accessors ((v v) (V V) (T1 T1) (H H) (f f) (F F) (g g) (G G) (Gp Gp) (logpgvgv logpgvgv)) filter
		(setf V (array-dimension v ROWS))
		(setf T1 (array-dimension v COLUMNS))
		(assert (and (> H 0) (> V 0) (> T1 0)))
		
		; The following five lines cannot be factored due to accessor restrictions
		(setf f (make-array T1 :initial-element (zeros H 1)))
		(setf F (make-array T1 :initial-element (zeros H H)))
		(setf g (make-array T1 :initial-element (zeros H 1)))
		(setf G (make-array T1 :initial-element (zeros H H)))
		(setf Gp (make-array (1- T1) :initial-element (zeros H H)))
		(setf logpgvgv (make-array T1 :initial-element 0))
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
; L (log likelihood of sequence log p(v(1:T1)))
(defmethod forward ((filter kalman-filter) &key (debug-p nil))
	(with-accessors ((T1 T1) (f f) (F F) (logpgvgv logpgvgv) (L L)) filter
		(when debug-p (print "Starting forward"))
		(loop for i below T1 do (forward-update filter i :debug-p debug-p))
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
	(with-accessors ((T1 T1) (f f) (F F) (g g) (G G) (Gp Gp)) filter
		(when debug-p (print "Starting backward"))
		(setf (aref g (1- T1)) (aref f (1- T1)))
		(setf (aref G (1- T1)) (aref F (1- T1)))
		(loop for i from (- T1 2) downto 0 do (backward-update filter i))
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
; Slots returned - modifications of g, G, Gp (smoothed mean, cov, cross moment <h_t h_{i+1}|v(1:T1)>)
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
	(with-accessors ((T1 T1) (v v) (f f) (F F) (g g) (G G) (Gp Gp)) filter
		;Refresh T1 in case synthetic observations added
		(setf T1 (array-dimension v COLUMNS)) 
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
			(V (array-dimension v ROWS)) (T1 (array-dimension v COLUMNS)))
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
						(loop for t1 below (- T1 1)
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
						(setf HH (%+% HH (aref G (1- T1)) (%*% (aref g (1- T1)) (tr (aref g (1- T1))))))
						(setf VH (%+% VH (%*% 
											(%-% (2d-aref-range v nil (1- T1)) meanV)
											(tr (aref g (1- T1)))))
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
						(loop for t1 below T1
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
							(setf meanV (%.*% (/ 1 T1) meanVtmp))
						)
						(setf CovV (%.*% (/ 1 T1) CovVtmp))
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
						(loop for t1 below (- T1 1)
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
						(setf CovH (%.*% CovHtmp (/ 1 (1- T1))))
						(setf CovH (%.*% .5 (%+% CovH (tr CovH))))
						(when diagCovH (setf CovH (diag CovH)))
						(when (not (loop for i below (array-total-size meanH) 
							always (zerop (aref meanH i 0)))) (setf meanH (%.*% (/ 1 (1- T1)) meanHtmp)))
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
	(let* ((T1 50)
			(lat (make-array `(1 ,T1) 
				:initial-contents (list (loop for i from 1 to T1 collect (sin (* .2 i))))))
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
								"original-lisp-data.csv"
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
								"reconstructed-lisp-data.csv"
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
	(let* ((T1 400) (Delta .1) (A (eye 6)) (B #2A((0d0 1d0 0d0 0d0 0d0 0d0) (0d0 0d0 0d0 1d0 0d0 0d0))))
		(loop for ind in (list '(0 4) '(1 0) '(2 5) '(3 2))
			do
			(setf (aref A (car ind) (cadr ind)) Delta)
		)
		(let ((h (zeros 6 T1)) (sigV 50d0) (sigH .00001) (v (zeros 2 T1))) ;FIXME Convert sigH to double
			(setf (aref h 1 0) (random 1d0))
			(setf (aref h 3 0) (random 1d0))
			(setf (aref h 0 0) (* 15 (random 1d0)))
			(setf (aref h 2 0) (* 15 (random 1d0)))
			(setf (aref h 4 0) (random 1d0))
			(setf (aref h 5 0) (- (random 1d0)))
			(2d-array-modify v (%+% (%*% B (2d-aref-range h nil 0)) (%.*% sigV (randmat 2 1))) nil 0)
			(loop for t1 from 1 to (1- T1)
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
								(format nil "~a.csv" mat-name))
					)
				)
			)
		)
	)
)

(defun demo-LDS-prediction (obs-num num-new-preds)
	(let* ((T1 obs-num)
			(lat (make-array `(1 ,T1) 
				:initial-contents (list (loop for i from 1 to T1 collect (sin (* .2 i))))))
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
						for filename in (list "noisy-lisp-data" "filtered-lisp-data")
						do
						(with-open-file (str (format nil 
												"~a.csv"
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