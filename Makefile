all: tsmc

tsmc: tsmc.c
	gcc tsmc.c -lm -o tsmc
