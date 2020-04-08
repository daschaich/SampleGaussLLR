# Quotation marks in shell commands below needed to handle "(" and ")"
CC = g++
CFLAGS = -O3 -Wall -std=c++11
LIBRARY_FLAGS = -lm
EXECUTABLES = 3d_U1_ora 4d_U1_mono

main: $(EXECUTABLES)

clean:
	rm -f $(EXECUTABLES) *.o

3d_U1_ora: metropolisu(1)overrelax.o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o 3d_U1_ora "metropolisu(1)overrelax.o"

metropolisu(1)overrelax.o:
	$(CC) $(CFLAGS) -c "metropolisu(1)overrelax.cpp"

4d_U1_mono: metropolisU(1)4dmonopoletest.o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o 4d_U1_mono "metropolisU(1)4dmonopoletest.o"

metropolisU(1)4dmonopoletest.o:
	$(CC) $(CFLAGS) -c "metropolisU(1)4dmonopoletest.cpp"
