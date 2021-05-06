all: functions.o main.o program
functions.o: functions.asm
	nasm -f elf32 -o functions.o functions.asm
main.o: main.c
	gcc -m32 -c -lm -o main.o main.c
program: main.o functions.o
	gcc -m32 -o program main.o functions.o
clean:
	rm -f *.o
	rm -f program
