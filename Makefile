# Configuration of the executable
TARGET = sph

INSTALL_PATH = $(PWD)

# Compiler configuration
CXX      = g++
CXXFLAGS = -Wall -Werror -Wextra -Wshadow -g
#CXXFLAGS = -Wall -Werror -Wextra -Wshadow -O3
#lFLAGS = -lOpenCL
LFLAGS = -L/home/adithya/Downloads/AMD-APP-SDK-v2.4-lnx64/lib/x86_64
IFLAGS = -I/boost_1_46_1/
COMP = -c 
INC = 
INC_LIB = -lOpenCL


OBJ = sph.o \
      param.o \
      particle.o \
    
default: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LFALGS) $(IFLAGS) -o $(TARGET) $(OBJ) $(INC) $(INC_LIB)

# clean
clean:
	@echo "clean ..."
	@$(RM) sph *.o *.md *.vtk
