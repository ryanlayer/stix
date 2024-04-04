BIN=../bin
OBJ=../obj
GIGGLE_ROOT=../../giggle
HTS_ROOT=$(GIGGLE_ROOT)/lib/htslib
# CFLAGS=-O2 -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized  
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
INCLUDES=-I$(HTS_ROOT) -I$(GIGGLE_ROOT) -I$(SQLITE_ROOT)
LIBS=-ldl -lz -lm -pthread -lcurl
LIBS+=-lcrypto
SQLITE_ROOT=../../sqlite-amalgamation-3170000
CC=gcc

LIB=ped.o \
    search.o \
	sharding_utils.o

.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))

.SECONDARY: $(OBJS)

# for debugging generates stix_dev
PROG_DEV=stix_dev
LIST_DEV=$(SQLITE_ROOT)/sqlite3.o $(addprefix $(BIN)/, $(PROG_DEV))
alldev: check-env $(LIST_DEV)


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

# for debugging generates stix_dev
$(BIN)/$(PROG): stix.c $(OBJS)
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


testshard:all
	../bin/stix -B ../test2/shard.tsv	-s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL

test:all
	../bin/stix -B ../test2/shard.tsv	-s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL


testvcf:all
	valgrind --tool=memcheck --leak-check=full  ../bin/stix -B ../test2/shard.tsv	-f ../test2/demo.vcf	-s 10 -v Sample

testvcf2:all
	../bin/stix -B ../test2/shard.tsv	-f ../test2/demo.vcf	-s 10 -v Sample


testQ:all
	../bin/stix -B ../test2/shard.tsv -Q ../test2/table.txt	-s 10

testj:all
	../bin/stix -j -B ../test2/shard.tsv	-s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL

testD:all
	../bin/stix -D -B ../test2/shard.tsv	-s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL

testa:all
	../bin/stix -a population,Sex,Alt_File -B ../test2/shard.tsv	-s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL -j

testsingle:all
	../bin/stix -a population,Sex,Alt_File   -i  ../../test/benchmark2/08.1.healthy_populations/03.1.giggle_idx_00/ -d   ../../test/benchmark2/08.1.healthy_populations/03.1.stix_idx_00.db      -s 10 -l 1:92345-92467 -r 1:116789-116799 -t DEL

