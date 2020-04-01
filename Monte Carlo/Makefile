CC = g++
CFLAGS = -O3 -Wall -std=c++11
INCLUDE_FLAGS =
LIBRARY_FLAGS = -lm
EXECUTABLE = U1_ora

OBJECT_FILES = "metropolisu(1)overrelax.o"

main: metropolisu(1)overrelax.o
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) $(LIBRARY_FLAGS) -o $(EXECUTABLE) $(OBJECT_FILES)

clean:
	rm -f $(EXECUTABLE) *.o

metropolisu(1)overrelax.o: metropolisu(1)overrelax.cpp
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c "metropolisu(1)overrelax.cpp"
