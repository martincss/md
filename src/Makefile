.PHONY: clean md

CC = gcc
LD = $(CC)
CFLAGS = -o3
LDFLAGS = -lm

SOURCE_C = $(wildcard *.c)
OBJECTS_C = $(patsubst %.c, %.o, $(SOURCE_C))
EXECUTABLE = ../bin/md.e

md: $(EXECUTABLE)

dynamics.o: dynamics.c
	$(CC) $(CFLAGS) -c $^ -o $@

../bin/md.e: dynamics.o md.c
	$(LD) $^ $(LD_FLAGS) -o $@

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE)
