CC = gcc
CFLAGS = -g -Wall -O3 -lm -std=c99 -fPIC
LD = gcc
LDFLAGS = -lm 
RM = /bin/rm -f
OBJS = stok1.o complex.o nrutil.o array.o mie.o
PROG = iquv

# top-level rule, to compile everything.
all: $(PROG)

# rule to link the program
$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -lm -o $(PROG) 

# rule for file "stok1.o".
stok1R.o: stok1R.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c stok1.c -lm -fPIC

# rule for file "nrutil.o".
nrutil.o: nrutil.c
	$(CC) $(CFLAGS) -c nrutil.c -lm -fPIC


# rule for file "complex.o".
complex.o: complex.c
	$(CC) $(CFLAGS) -c complex.c -lm -fPIC


# rule for file "mie.o".
mie.o: mie.c
	$(CC) $(CFLAGS) -c mie.c -lm -fPIC

	
# rule for file "array.o".
array.o: array.c
	$(CC) $(CFLAGS) -c array.c -lm -fPIC



# rule for cleaning re-compilable files.
clean:
	$(RM) $(PROG) $(OBJS)
