.PHONY: clean md

# variables (todo lo rojo)
CC = gcc
LD = $(CC)
CFLAGS = -std=gnu99 -O3
LD_FLAGS = -lm

SOURCE_C = $(wildcard *.c)					# todas las cosas que terminan en '.c'
OBJECTS_C = $(patsubst %.c, %.o, $(SOURCE_C))
EXECUTABLE = ../bin/md.e						# defino el ejecutable

md: $(EXECUTABLE)

dynamics.o: dynamics.c							# target: dependencia
	$(CC) $(CFLAGS) -c $^ -o $@			# reglas de compilacion

# $@ = target
# $^ = todas las dependencias

../bin/md.e: dynamics.o md.c				# archivo: dependencias
	$(LD) $^ $(LD_FLAGS) -o $@

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE)
