$(call path_prepend,$(NETCDF_DIR)/lib/pkgconfig,PKG_CONFIG_PATH)

NETCDF_DIR != $(PKG_CONFIG) --variable=prefix netcdf

ifeq ($(NETCDF_DIR),)
  $(call cwarning,Netcdf root directory NETCDF_DIR -- undefined)
else
  $(call cinfo,-- Netcdf root directory NETCDF_DIR -- $(NETCDF_DIR))
endif

ifeq ($(wildcard $(NETCDF_DIR)/include/netcdf_par.h),)
  ifeq ($(wildcard $(NETCDF_DIR)/include/netcdf.h),)
    $(call cwarning,Include file $$NETCDF_DIR/include/netcdf.h -- not found)
    undefine DEFINE_PARALLEL_NETCDF
  else
    $(call cinfo,-- Include file $$NETCDF_DIR/include/netcdf.h -- found)
    $(call cinfo,--     Using sequential Netcdf)
    DEFINE_PARALLEL_NETCDF = \#undef USE_PARALLEL_NETCDF
  endif
else
  $(call cinfo,-- Include file $$NETCDF_DIR/include/netcdf_par.h -- found)
  $(call cinfo,--     Using parallel Netcdf)
  DEFINE_PARALLEL_NETCDF = \#define USE_PARALLEL_NETCDF
endif
