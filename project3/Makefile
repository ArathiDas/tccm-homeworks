# Compiler and Flags
CC = gcc
CFLAGS = -Wall -g -O2

LIBS = -lm

# Executable and Object Files
TARGET = dynamics
OBJS = dynamics.o utils.o error.o


# Default Target
all: $(TARGET)

# Rule to Build the Executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Compile Source Files into Object Files
dynamics.o: dynamics.c utils.h error.h
	$(CC) $(CFLAGS) -c dynamics.c $(LIBS)

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c $(LIBS)

error.o: error.c error.h
	$(CC) $(CFLAGS) -c error.c $(LIBS)

# Clean Target to Remove Build Artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Run the Executable
run: all
	./$(TARGET)

