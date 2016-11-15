$(call path_prepend,$(NETCDF_DIR)/lib/pkgconfig,PKG_CONFIG_PATH)

NETCDF_DIR != $(PKG_CONFIG) --variable=prefix netcdf

ifeq ($(NETCDF_DIR),)
  $(call cwarning,     Pkg-config file netcdf.pc -- not found)
endif
ifeq ($(wildcard $(NETCDF_DIR)/include/netcdf_par.h),)
  define msg
     Include file $$(NETCDF_DIR)/include/netcdf_par.h -- not found.
       Install parallel Netcdf4/HDF5 and/or set NETCDF_DIR?
  endef
  $(call cwarning,$(msg))
endif
