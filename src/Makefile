CC=g++
CFLAGS= -O3 -Wall -std=c++11 -I/usr/local/Cellar/igraph/0.7.1_6/include/igraph 
LDFLAGS= -lm -Wall -L/usr/local/lib -ligraph -lglpk 
EXEC= spectral post_processing modularity optimal optimal_community

all: $(EXEC)

spectral : graph.o spectral.o 
	$(CC) -o $@ $^ $(LDFLAGS)
modularity : graph.o modularity.o
	$(CC) -o $@ $^ $(LDFLAGS)
post_processing : graph.o post_processing.o 
	$(CC) -o $@ $^ $(LDFLAGS)
optimal : graph.o optimal.o 
	$(CC) -o $@ $^ $(LDFLAGS) 	
optimal_community : graph.o optimal.o 
	$(CC) -o $@ $^ $(LDFLAGS) 	

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CC) -o $@ -c -g $< $(CFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c -g $< $(CFLAGS)

clean:
	rm -f *.o *~ $(EXEC)
