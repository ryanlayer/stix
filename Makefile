BIN=bin
OBJ=obj

all: 
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE)

debug:

	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE) -f Makefile-debug

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
