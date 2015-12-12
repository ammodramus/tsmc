all: tsmc

tsmc:
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c -lm -o tsmc

debug:
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c -lm -o tsmc -g
