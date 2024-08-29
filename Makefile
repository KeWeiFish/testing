# Define compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11 -Wall

# Name of the executable
TARGET1 = main
TARGET2 = test
.PHONY: all clean $(TARGET1) 
# Build rule
all: $(TARGET1) 

$(TARGET1): $(TARGET1).cpp
	$(CC) $(CFLAGS) $(TARGET1).cpp -o $(TARGET1)


clean:
	rm -f $(TARGET1) $(TARGET1).o
	rm -rf GBresults
	rm -rf *Convergence*