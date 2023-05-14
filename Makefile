all:
	@echo "This project implemented with CMake"
	@echo "Try this:"
	@echo mkdir build
	@echo cd build
	@echo cmake ..
	@echo make -j 8

prepare:
	rm -rf build
	mkdir build