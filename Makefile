all: tsmc

tsmc: tsmc.c
	gcc tsmc.c hmm.c definitions.c qtscases.c -lm -o tsmc

debug:
	gcc tsmc.c hmm.c definitions.c qtscases.c -lm -o tsmc -g
