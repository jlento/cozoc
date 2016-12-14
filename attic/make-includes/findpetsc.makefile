$(call path_prepend,$(PETSC_DIR)/lib/pkgconfig,PKG_CONFIG_PATH)

PETSC_DIR != $(PKG_CONFIG) --variable=prefix PETSc

ifeq ($(PETSC_DIR),)
  define msg
     Pkg-config file PETSc.pc -- not found
       Install PETSc and/or set variable PETSC_DIR?
  endef
  $(call cwarning,$(msg))
endif
