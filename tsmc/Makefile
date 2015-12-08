all: tsmc

tsmc:
	gcc tsmc.c hmm.c definitions.c qtscases.c -lm -o tsmc

debug:
	gcc tsmc.c hmm.c definitions.c qtscases.c -lm -o tsmc -g

simt3t2: simt3t2.c simt3t2.h random.c random.h
	gcc simt3t2.c random.c -O2 -lm -o simt3t2

