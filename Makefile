# Compiler
CC = mpicc

# Flags for the compiler
CFLAGS = 

# Flags for the linker
LDFLAGS = -lm

# Target executable
TARGET = mandelbrot_mw

# Source files
SRC = mandelbrot_mw.c

# Build the target
$(TARGET): $(SRC)
	$(CC) -o $(TARGET) $(SRC) $(CFLAGS) $(LDFLAGS)

# Clean up build files
clean:
	rm -f $(TARGET)
