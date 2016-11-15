(setq inhibit-splash-screen t
      initial-scratch-message nil)
(setq make-backup-files nil)
(defalias 'yes-or-no-p 'y-or-n-p)
(when window-system (set-frame-size (selected-frame) 100 24))
(when window-system (add-to-list 'default-frame-alist '(fullscreen . fullheight)))

(require 'cl)
(load "package")
(package-initialize)
(add-to-list 'package-archives
             '("melpa" . "http://melpa.milkbox.net/packages/") t)
(defvar jlento/packages '(magit
			  paredit
			  ggtags
			  sr-speedbar
			  auto-complete)
  "Default packages")

(defun jlento/packages-installed-p ()
  (loop for pkg in jlento/packages
        when (not (package-installed-p pkg)) do (return nil)
        finally (return t)))

(unless (jlento/packages-installed-p)
  (message "%s" "Refreshing package database...")
  (package-refresh-contents)
  (dolist (pkg jlento/packages)
    (when (not (package-installed-p pkg))
      (package-install pkg))))


;; (global-semantic-idle-summary-mode 1)

(require 'helm-config)
(helm-mode 1)

(require 'auto-complete-config)
(ac-config-default)

(require 'ggtags)
(add-hook 'c-mode-common-hook
          (lambda ()
            (when (derived-mode-p 'c-mode 'c++-mode 'java-mode)
              (ggtags-mode 1))))

(require 'sr-speedbar)
(setq
 sr-speedbar-right-side nil
 sr-speedbar-width 20
 sr-speedbar-width-console 20
 sr-speedbar-max-width 10
 sr-speedbar-delete-windows t)
(sr-speedbar-open)

(require 'whitespace)
(setq-default show-trailing-whitespace t)

;; (add-hook 'c-mode-hook 
;; 	  '(lambda () 
;; 	     (gtags-mode t)
;; 	     ))
;; (defun gtags-root-dir ()
;;   "Returns GTAGS root directory or nil if doesn't exist."
;;   (with-temp-buffer
;;     (if (zerop (call-process "global" nil t nil "-pr"))
;; 	(buffer-substring (point-min) (1- (point-max)))
;;       nil)))
;; (defun gtags-update ()
;;   "Make GTAGS incremental update"
;;   (call-process "global" nil nil nil "-u"))
;; (defun gtags-update-hook ()
;;   (when (gtags-root-dir)
;;     (gtags-update)))
;; (add-hook 'after-save-hook #'gtags-update-hook)



;;(when window-system 
;;  (speedbar t))
