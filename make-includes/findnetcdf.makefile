$(call path_prepend,$(NETCDF_DIR)/lib/pkgconfig,PKG_CONFIG_PATH)

NETCDF_DIR != $(PKG_CONFIG) --variable=prefix netcdf

ifeq ($(NETCDF_DIR),)
  $(call cwarning,Netcdf root directory NETCDF_DIR -- undefined)
else
  ifndef SILENT
    $(info -- Netcdf root directory NETCDF_DIR -- $(NETCDF_DIR))
  endif
endif

ifeq ($(wildcard $(NETCDF_DIR)/include/netcdf_par.h),)
  undefine USE_PARALLEL_NETCDF
  ifeq ($(wildcard $(NETCDF_DIR)/include/netcdf.h),)
    $(call cwarning,Include file $$NETCDF_DIR/include/netcdf.h -- not found)
  else
    $(call cinfo,-- Include file $$NETCDF_DIR/include/netcdf.h -- found)
    $(call cinfo,--     Using sequential Netcdf)
  endif
else
  USE_PARALLEL_NETCDF = 1
  $(call cinfo,-- Include file $$NETCDF_DIR/include/netcdf_par.h -- found)
  $(call cinfo,--     Using parallel Netcdf)
endif
