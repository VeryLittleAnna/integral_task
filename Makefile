M1 = -DCOMBINE
M2 = 
KEY = -D1
ifeq ($(METHOD), COMBINE)
	KEY = $(M1)
else
	KEY = $(M2)
endif
all: functions.o main.o program
functions.o: functions.asm
	nasm -f elf32 -o functions.o functions.asm
main.o: main.c
	gcc -m32 -c -lm $(KEY) -o main.o main.c
program: main.o functions.o
	gcc -m32 -o program main.o functions.o
clean:
	rm -f *.o
	rm -f program
