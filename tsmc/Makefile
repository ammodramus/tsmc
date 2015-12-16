all: optimize

optimize:
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c em.c asa047.c -lm -o tsmc -O3 -DNDEBUG -march=native

debug:
	gcc tsmc.c hmm.c definitions.c qtscases.c data.c em.c asa047.c -lm -o tsmc -g

emissionstest:
	gcc emissionstest.c hmm.c definitions.c qtscases.c data.c em.c asa047.c -lm -o emissionstest -g

