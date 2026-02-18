TARGET = run

CXX = g++

CXXFLAGS = -Wall -Wextra -std=c++17 -g -fsanitize=address -O0

SRC = main.cpp \
      Cell/Cell.cpp \
      Phase/Phase.cpp \
      Components/Components.cpp \
      Face/Face.cpp \
      Mesh/Mesh.cpp \
      MUSCL/MUSCL.cpp \
      RiemanSolver/RiemanSolver.cpp \
      HyperbolicOperator/HyperbolicOperator.cpp \
      RelaxationOperator/RelaxationOperator.cpp \
      Solver/Solver.cpp

OBJ = $(SRC:.cpp=.o)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJ) $(TARGET)