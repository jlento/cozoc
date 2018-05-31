#!/bin/bash

# Bash utility functions

require () {
    local f rcode=0
    for f in "$@"; do
        F="$(tr '.' '_' <<<${f^^})"
        test -z "${!F}" && printf -v "$F" "%s" "$f"
        which ${!F} > /dev/null || {
            echo "Executable '$f' not found." >&2
            echo "    Try 'export $F=/full/path/to/$f'" >&2
            rcode=1
        }
    done
    return $rcode
}
