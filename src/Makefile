BIN=../bin
OBJ=../obj
GIGGLE_ROOT=../../giggle
HTS_ROOT=$(GIGGLE_ROOT)/lib/htslib
CFLAGS=-O2 -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized 
#CFLAGS=-g -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
INCLUDES=-I$(HTS_ROOT) -I$(GIGGLE_ROOT) -I$(SQLITE_ROOT)
LIBS=-ldl -lz -lm -pthread -lcurl
LIBS+=-lcrypto
SQLITE_ROOT=../../sqlite-amalgamation-3170000
CC=gcc

LIB=ped.o \
    search.o

.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))

.SECONDARY: $(OBJS)


PROG=stix

LIST=$(SQLITE_ROOT)/sqlite3.o $(addprefix $(BIN)/, $(PROG))

all: check-env $(LIST)

$(SQLITE_ROOT)/sqlite3.o: $(SQLITE_ROOT)/sqlite3.c
	$(CC) -c $(SQLITE_ROOT)/sqlite3.c -o $(SQLITE_ROOT)/sqlite3.o

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< \
		$(INCLUDES) 

$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		$(INCLUDES) \
		-DSAMTOOLS=1 \
		$(GIGGLE_ROOT)/lib/libgiggle.a \
		$(HTS_ROOT)/libhts.a \
		$(SQLITE_ROOT)/sqlite3.o \
		$(LIBS)

check-env:
ifdef TRAVIS_BUILD_DIR
HTS_ROOT=$(TRAVIS_BUILD_DIR)/htslib
GIGGLE_ROOT=$(TRAVIS_BUILD_DIR)/giggle
endif
