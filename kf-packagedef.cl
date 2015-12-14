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
;; $Revision: 1 $
;;
;; -=End Copyright Notice=-

(in-package :cl-user)

(provide "kf")

(eval-when (load compile eval)
  ;(require :2is-build-support "S:\\_Vault_Working\\Libraries\\2isBuild\\src\\development\\trunk\\build.cl")
  ;(pcheck-module "statistics" "s:\\_vault_working\\libraries\\statistics\\src\\development\\trunk\\statistics.lpr")
  ;(pcheck-module "acLib" "S:\\_Vault_Working\\Libraries\\allegrocachelib\\branches\\1.7\\aclib.lpr")
)

(defpackage :statistics)
(defpackage :db.allegrocache)

(defpackage :kf
  (:use :cl :excl :db.allegrocache)
  (:export))
   