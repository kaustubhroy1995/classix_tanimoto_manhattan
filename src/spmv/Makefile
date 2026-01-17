# Name of the file to be compiled
FILENAME := spsubmatxvec

# Required C and CXX compiler
CXX := clang
CC := clang

# Target Rules
all:build

build:srcfile

srcfile:
	$(CC) -fPIC -O3 -o $(FILENAME).o -c $(FILENAME).c
	$(CC) --shared -fPIC $(FILENAME).o -o $(FILENAME).so


run: build
	$(EXEC) ./$(FILENAME)

clean:
	rm -rf *o
	rm -rf $(FILENAME).so