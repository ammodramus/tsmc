all: optimize

optimize: *.[ch]
	gcc tsmc.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o tsmc -O3 -DNDEBUG -march=native

debug: *.[ch]
	gcc tsmc.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o tsmc_debug -g

debugmsg: *.[ch] 
	gcc tsmc.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o tsmc -O3 -march=native

profile: *.[ch]
	gcc tsmc.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o tsmcprofile -pg

emissionstest: *.[ch]
	gcc emissionstest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o emissionstest -g

inputtest: *.[ch]
	gcc inputtest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o inputtest -g

forwardtest: *.[ch]
	gcc forwardtest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o forwardtest -g

backwardtest: *.[ch]
	gcc backwardtest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o backwardtest -g

normtest: *.[ch]
	gcc normtest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o normtest -g

pitest: *.[ch]
	gcc pitest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o pitest -g

gammatest: *.[ch]
	gcc gammatest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o gammatest -g

expecttest: *.[ch]
	gcc expecttest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o expecttest -g

logliketest: *.[ch]
	gcc likelihoodtest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c util.c -lm -o logliketest -g

qtstest: *.[ch]
	gcc qtstest.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c options.c util.c -lm -o qtstest -g

optionstest: *.[ch]
	gcc optionstest.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c -lm -o optionstest -g

optimtest: *.[ch]
	gcc optimtest.c definitions.c asa047.c util.c -lm -o optimtest -g

diptriptestpts: *.[ch]
	gcc diptriptestpts.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o diptriptestpts -g

diptriptestsurface: *.[ch]
	gcc diptriptestsurface.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o diptriptestsurface -O3 -DNDEBUG -march=native

memprofile: *.[ch]
	gcc tsmc.c options.c hmm.c definitions.c qtscases.c qtssuppdt.c data.c em.c asa047.c random.c util.c -lm -o tsmc -O3 -DDEBUGMEMORY -march=native

testgenlikes: *.[ch]
	gcc testgenlikedata.c genotypelikelihoods.c definitions.c util.c -lm -g -o testgenlikedata
