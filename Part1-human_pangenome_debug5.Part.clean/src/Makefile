.PHONY: release clean debug
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
THREADS := 4
export LIBCUCKOO = -I${ROOT_DIR}/lib/libcuckoo
.DEFAULT_GOAL := all
CXX :=  g++
CXXFLAGS += -Wall -Wextra -pthread -std=c++17 -g
CXXFLAGS += -Wno-missing-field-initializers
CXXFLAGS += -Ilibcuckoo/
BIN_OUT := ${BIN_DIR}/pangenomeDebug5.Part.clean

LDFLAGS += -pthread -std=c++17 -rdynamic -lz -lcurl

debug: CXXFLAGS += -g
debug: PANGENOME

release: PANGENOME



mainObj := ${patsubst %.cpp, %.o, ${wildcard *.cpp}}

%.o: %.cpp parallel.h progress_bar.h logger.h config.h
	${CXX} -c ${CXXFLAGS} $< -o $@


PANGENOME: ${mainObj}
	${CXX} ${mainObj} -o ${BIN_OUT} ${LDFLAGS}


clean:
	rm -f ${mainObj}
	rm -f ${BIN_OUT}
