ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
THREADS := 4

.PHONY: all debug clean
.DEFAULT_GOAL := all

export BIN_DIR = ${ROOT_DIR}/bin


all: 
	make release -C src -j ${THREADS}
	# cp ${BIN_DIR}/pangenomeDebug5.Part.clean /work/path
debug:
	make debug -C src -j ${THREADS}
clean: 
	make clean -C src
	rm -f ${BIN_DIR}/pangenomeDebug5.Part
	rm /work/path
