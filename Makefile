# Quotation marks in shell commands below needed to handle "(" and ")"
CC = g++
CFLAGS = -O3 -Wall -std=c++11
LIBRARY_FLAGS = -lm
EXECUTABLES = 3d_U1_LLR

main: $(EXECUTABLES)

clean:
	rm -f $(EXECUTABLES) *.o

3d_U1_LLR: LLRU(1).o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o 3d_U1_ora "LLRU(1).o"

LLRU(1).o:
	$(CC) $(CFLAGS) -c "LLRU(1).cpp"
