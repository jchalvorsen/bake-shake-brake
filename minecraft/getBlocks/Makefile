CC = g++
PROG = getBlocks

LIBS	= -lz

CFLAGS = -Wall -g 

SRCS = getBlocks.cpp zlib.cpp 


all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)
