CFLAGS=gcc -std=c99 -O3 -Wall
LDFLAGS=-lm

all: linked_list energy cells code
	$(CFLAGS) build/linked_list.o build/energy.o build/cells.o build/code.o -o build/start $(LDFLAGS)

cells:
	$(CFLAGS) -c src/cells.c -o build/cells.o $(LDFLAGS)

energy:
	$(CFLAGS) -c src/energy.c -o build/energy.o $(LDFLAGS)

linked_list:
	$(CFLAGS) -c src/linked_list.c -o build/linked_list.o $(LDFLAGS)

code:
	$(CFLAGS) -c src/code.c -o build/code.o $(LDFLAGS)
