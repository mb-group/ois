# Declaration of variables
CC = mpic++
CFLAGS = -O3 -mtune=native  -std=c++17

# File names
EXEC  = orthoseq
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
INCLUDE = -Iinclude/ 

# Main target
default:all

all:  	$(OBJECTS)
	$(CC) ${INCLUDE} $(OBJECTS)  -O3  -std=c++17 -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) $(CFLAGS) ${INCLUDE}  -c   $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
	clear
