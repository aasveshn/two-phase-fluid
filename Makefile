TARGET = run

CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

SRC = main.cpp  Cell/Cell.cpp 
OBJ = $(SRC:.cpp =.o)

$(TARGET): $(OBJ)
	$(CXX) $(CXXfLAGS) -o $@ $^
%.0: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean: 
	rm -f $(OBJ) $(TARGET)
