CPPCOMP = g++
DEBUGGER = gdb
OPT = -O3
#OPT = -g

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

program.exe: $(SRC_DIR)/*.cpp
	$(CPPCOMP) $(OPT) -I $(INCLUDE_DIR) $(SRC_DIR)/*.cpp -o ${BUILD_DIR}/program.exe

run: program.exe
	./${BUILD_DIR}/program.exe

clean:
	rm -f ${BUILD_DIR}/*.exe

debug:
	$(DEBUGGER) ./${BUILD_DIR}/program.exe