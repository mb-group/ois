# Declaration of variables
CC = mpic++
CFLAGS = -O3 -mtune=native -march=native -std=c++17

# File names
EXEC  = c3d
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
INCLUDE = -Iinclude/ 

# Main target
default:all

all:  	$(OBJECTS)
	$(CC) ${INCLUDE} $(OBJECTS)  -O3  -std=c++11 -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) $(CFLAGS) ${INCLUDE}  -c   $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
	clear
