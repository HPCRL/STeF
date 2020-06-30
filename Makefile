EXE_DIR = exe
BIN_DIR = bin
SRC_DIR = src
OBJ_DIR = obj
DIRS = $(OBJ_DIR) $(BIN_DIR)

CC = icpc

EXE = $(wildcard $(EXE_DIR)/*.c)
BIN = $(EXE:$(EXE_DIR)/%.c=$(BIN_DIR)/%.exe)
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CPPFLAGS += -I /home/emre/papi/src/ -I /home/semrekurt/Desktop/code/blis/include/haswell/ -D_POSIX_C_SOURCE=200112L -DBLIS_VERSION_STRING="0.7.0-10"
CFLAGS += -Wall -Wno-write-strings -g -std=c++11 -O3 -Ofast  -qopenmp $(EXTRA) -march=native -restrict -mkl
#CFLAGS += -Wall -Wno-write-strings -g -std=c++11 -O0 -qopenmp $(EXTRA) -march=native -restrict -mkl
LDFLAGS += -L /home/aravind/bin/papi/lib -L home/semrekurt/Desktop/code/blis/lib/haswell/libblis.a -lm -lpthread
LDLIBS += #-lpapi


.PHONY: all clean

all: | ${DIRS} $(BIN)

papi: EXTRA += "-D PAPI"
papi: clean dir 
papi: $(BIN)

${DIRS}:
	mkdir $@

$(BIN_DIR)/%.exe: $(OBJ) $(EXE_DIR)/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $^  -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	$(RM) -rf $(DIRS)
