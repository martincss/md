.PHONY: clean md 1a 2 3

# variables (todo lo rojo)
CC = gcc
LD = $(CC)
CFLAGS = -std=c99 -O3  # -std=c99
LD_FLAGS = -lm

SOURCE_C = $(wildcard *.c)					# todas las cosas que terminan en '.c'
OBJECTS_C = $(patsubst %.c, %.o, $(SOURCE_C))
EXECUTABLE = ../bin/md.exe						# defino el ejecutable
1a = ../bin/md_1a.exe
2 = ../bin/md_2.exe
3 = ../bin/md_3.exe

md: $(EXECUTABLE)

1a: $(1a)

2: $(2)

3: $(3)

dynamics.o: dynamics.c							# target: dependencia
	$(CC) $(CFLAGS) -c $^ -o $@			# reglas de compilacion

pruebas_pos.o: pruebas_pos.c							# target: dependencia
	$(CC) $(CFLAGS) -c $^ -o $@			# reglas de compilacion

# $@ = target
# $^ = todas las dependencias

../bin/md.exe: dynamics.o pruebas_pos.o md.c				# archivo: dependencias
	$(LD) $^ $(LD_FLAGS) -o $@

../bin/md_1a.exe: dynamics.o pruebas_pos.o md_1a.c				
	$(LD) $^ $(LD_FLAGS) -o $@

../bin/md_2.exe: dynamics.o pruebas_pos.o md_2.c				
	$(LD) $^ $(LD_FLAGS) -o $@	

../bin/md_3.exe: dynamics.o pruebas_pos.o md_3.c				
	$(LD) $^ $(LD_FLAGS) -o $@	

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE) $(1a)
