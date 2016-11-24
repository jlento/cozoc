# Constants

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
define LF


endef

RED    != tput setaf 1
GREEN  != tput setaf 2
YELLOW != tput setaf 3
RESET  != tput sgr0

# Functions

cerror = $(error $(RED)$1$(RESET))

ifeq ($(filter s,$(MAKEFLAGS)),)
    cwarning = $(warning $(YELLOW)$1$(RESET))
    cinfo    = $(info $(GREEN)$1$(RESET))
else
    cwarning =
    cinfo =
endif

_path_prepend = $(eval $2 := $(subst $(SPACE),:,$(strip $1 $($2))))
path_prepend  = $(if $(wildcard $1),$(call _path_prepend,$1,$2))


# Project defaults

$(call path_prepend,$(CURDIR)/lib/pkgconfig,PKG_CONFIG_PATH)


#.PHONY : all
# all : ; echo $(PKG_CONFIG_PATH)
