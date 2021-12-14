EXE_DIR = exe
BIN_DIR = bin
SRC_DIR = src
OBJ_DIR = obj
INC_DIR = inc

DIRS = $(OBJ_DIR) $(BIN_DIR)

CC = g++

LIKWID_DIR=/opt/likwid/
LIKWID_INC=$(LIKWID_DIR)/include
LIKWID_LIB=$(LIKWID_DIR)/lib

EXE = $(wildcard $(EXE_DIR)/*.cpp)
BIN = $(EXE:$(EXE_DIR)/%.cpp=$(BIN_DIR)/%.exe)
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CPPFLAGS += -I/usr/local/include   
CFLAGS = -Wall -Wno-write-strings -g -std=c++14 -O3 $(EXTRA) -march=native  -fopenmp -funroll-loops -fstrict-aliasing -fPIC -I/$(LIKWID_INC)
#CFLAGS += -Wall -Wno-write-strings -g -std=c++11 -O0 -fopenmp $(EXTRA) 
LDFLAGS += -L/usr/local/lib -L/$(LIKWID_LIB) -lpthread 
LDLIBS += -llikwid
#CFLAGS += -D OMP
CPPFLAGS += -D LIKWID_PERFMON

.PHONY: all clean

all: | ${DIRS} $(OBJ) $(BIN)

papi: EXTRA += "-D PAPI"
papi: clean dir 
papi: $(BIN)

${DIRS}:
	mkdir $@

$(BIN_DIR)/%.exe: $(OBJ) $(EXE_DIR)/%.cpp
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $^  -o $@ $(LDLIBS)

$(OBJ_DIR)/%.o:  $(SRC_DIR)/%.cpp 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ $(LDLIBS)

clean: 
	$(RM) -rf $(DIRS)
