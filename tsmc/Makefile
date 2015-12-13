all: tsmc

tsmc: *.[ch]
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c em.c -lm -o tsmc -O2 -DNDEBUG

debug:
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c em.c -lm -o tsmc -g
