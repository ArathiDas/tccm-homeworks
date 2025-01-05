# Compiler and Flags
CC = gcc
CFLAGS = -Wall -g -O2

LIBS = -lm

# Executable and Object Files
TARGET = dynamics
OBJS = src/dynamics.o src/utils.o src/error.o


# Default Target
all: $(TARGET)

# Rule to Build the Executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Compile Source Files into Object Files
dynamics.o: dynamics.c utils.h error.h
	$(CC) $(CFLAGS) -c src/dynamics.c $(LIBS)

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c src/utils.c $(LIBS)

error.o: error.c error.h
	$(CC) $(CFLAGS) -c src/error.c $(LIBS)

# Clean Target to Remove Build Artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Run the Executable
run: all
	./$(TARGET)

