;;; init.el -- EMACS configuration for C development

;;; Commentary:
;;;
;;; Excellent Youtube videos about using EMACS in C development:
;;;   https://www.youtube.com/watch?v=HTUE03LnaXA
;;; inspired this setup, too!

;;; Code:

;; Basic preferences
(custom-set-variables
      '(inhibit-splash-screen t)
      '(initial-scratch-message nil)
      '(indent-tabs-mode nil)
      '(c-basic-offset 4)
      '(make-backup-files nil))
(defalias 'yes-or-no-p 'y-or-n-p)
(column-number-mode)
(global-set-key (kbd "C-c a") 'align-current)
(when window-system (set-frame-width (selected-frame) 120))

;;; Boot use-package
(require 'package)
(setq package-enable-at-startup nil)
(add-to-list 'package-archives
             '("melpa" . "http://melpa.org/packages/"))
(package-initialize)
(unless (package-installed-p 'use-package)
  (package-refresh-contents)
  (package-install 'use-package))
(eval-when-compile
  (require 'use-package))
(require 'diminish)
(require 'bind-key)

(setq use-package-verbose t)


;;; The packages...

(use-package sr-speedbar
  :ensure t
  :config
  (custom-set-variables
   '(sr-speedbar-right-side nil)
   '(sr-speedbar-width 20)
   '(sr-speedbar-width-console 20)
   '(sr-speedbar-max-width 10)
   '(sr-speedbar-delete-windows t))
  (sr-speedbar-open))

(use-package auto-complete
  :ensure t
  :config
  (require 'auto-complete-config)
  (add-to-list 'ac-dictionary-directories
               (expand-file-name
                "~/.emacs.d/elpa/auto-complete-20161029.643/dict"))
  (setq ac-comphist-file
        (expand-file-name "~/.emacs.d/ac-comphist.dat"))
  (ac-config-default))

(use-package auto-complete-c-headers
  :ensure t
  :config
  (defun my:ac-c-header-init ()
    (add-to-list 'ac-sources 'ac-source-c-headers)
    (add-to-list 'achead:include-directories '"/usr/include")
    (add-to-list 'achead:include-directories
                 '"/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real/include"))
  (add-hook 'c-mode-hook 'my:ac-c-header-init))

(use-package yasnippet
  :ensure t
  :config
  (yas-global-mode 1))

(use-package iedit
  :ensure t
  :config
  (define-key global-map (kbd "C-c ;") 'iedit-mode))

(use-package semantic
  :ensure t
  :config
  (global-semantic-decoration-mode t)
  (global-semantic-idle-scheduler-mode t)
  (global-semantic-idle-completions-mode t)
  (global-semantic-idle-summary-mode t)
  (semantic-mode t)
  (require 'semantic/ia)
  (require 'semantic/bovine/gcc)
  (defun my:add-semantic-to-autocomplete()
    (add-to-list 'ac-sources 'ac-source-semantic))
  (add-hook 'c-mode-hook 'my:add-semantic-to-autocomplete))

(message "JEEEEEEEEEEEEEE")
(message (if (getenv "PETSC_DIR")
             (concat (getenv "PETSC_DIR") "/include")
           (shell-command-to-string "pkg-config --variable=includedir PETSc")))
(message (car (split-string (shell-command-to-string
                             "pkg-config --variable=includedir mpi"))))
(message "JEEEEEEEEEEEEEE")

(use-package ede
  :ensure t
  :config
  (global-ede-mode t)
  (ede-cpp-root-project
   "cozoc"
   :name "cozoc"
   :file load-file-name
   :include-path '("/src")
   :system-include-path
   (list "/usr/include"
         (car (split-string (shell-command-to-string
                             "pkg-config --variable=includedir mpi")))
         (if (getenv "NETCDF_DIR")
             (concat (getenv "NETCDF_DIR") "/include")
           (car (split-string (shell-command-to-string
                               "pkg-config --variable=includedir netcdf"))))
         (if (getenv "PETSC_DIR")
             (concat (getenv "PETSC_DIR") "/include")
           (car (split-string (shell-command-to-string
                               "pkg-config --variable=includedir PETSc")))))))

(use-package whitespace
  :init
  (dolist (hook '(prog-mode-hook text-mode-hook))
    (add-hook hook #'whitespace-mode))
  (add-hook 'before-save-hook #'whitespace-cleanup)
  :config
  (setq whitespace-line-column 80) ;; limit line length
  (setq whitespace-style '(face tabs empty trailing lines-tail)))

(defun c-reformat-buffer ()
  (interactive)
  (shell-command-on-region
   (point-min)
   (point-max)
   "bash -c 'clang-format $1 | astyle -d --style=lisp' --"
   (buffer-name) t))
(define-key global-map (kbd "C-c r") 'c-reformat-buffer)


;; Fix iedit bug in Mac

;; Flycheck

;; start flymake-google-cpplint-load
;; let's define a function for flymake initialization
;; (defun my:flymake-google-init ()
;;   (require 'flymake-google-cpplint)
;;   (custom-set-variables
;;    '(flymake-google-cpplint-command "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/cpplint"))
;;   (flymake-google-cpplint-load)
;; )
;; (add-hook 'c-mode-hook 'my:flymake-google-init)
;; (add-hook 'c++-mode-hook 'my:flymake-google-init)

;; ; start google-c-style with emacs
;; (require 'google-c-style)
;; (add-hook 'c-mode-common-hook 'google-set-c-style)
;; (add-hook 'c-mode-common-hook 'google-make-newline-indent)


;; Semantic
;; (semantic-mode 1)
;; (defun my:add-semantic-to-autocomplete()
;;   (add-to-list 'ac-sources 'ac-source-semantic))
;; (add-hook 'c-mode-common-hook 'my:add-semantic-to-autocomplete)

;; ;; ede mode
;; (global-ede-mode 1)
;; (ede-cpp-root-project "my project" :file "src/cozoc.c"
;; 		      :include-path '("src" "/usr/include"))
;; ; you can use system-include-path for setting up the system header file locations.
;; ; turn on automatic reparsing of open buffers in semantic
;; (global-semantic-idle-scheduler-mode 1)
