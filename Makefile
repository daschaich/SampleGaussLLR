# Quotation marks in shell commands below needed to handle "(" and ")"
CC = g++
CFLAGS = -O3 -Wall -std=c++11
LIBRARY_FLAGS = -lm
EXECUTABLES = 3d_U1_LLR 4d_U1_LLR gauss_LLR

main: $(EXECUTABLES)

clean:
	rm -f $(EXECUTABLES) *.o

3d_U1_LLR: LLRU(1)newtry.o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o 3d_U1_LLR "LLRU(1)newtry.o"

LLRU(1)newtry.o:
	$(CC) $(CFLAGS) -c "LLRU(1)newtry.cpp"

4d_U1_LLR: LLR.o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o 4d_U1_LLR LLR.o

LLR.o:
	$(CC) $(CFLAGS) -c LLR.cpp

gauss_LLR: LLRgaussc++.o
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o gauss_LLR "LLRgaussc++.o"

LLRgaussc++.o:
	$(CC) $(CFLAGS) -c "LLRgaussc++.cpp"
