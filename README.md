# How to compile and run

```
mkdir build
cd build
cmake ..
make 
./bin/main
```

# How to compile and run on Windows

Download and install MSYS2 (https://www.msys2.org/) using default settings. Start the MSYS2 MINGW 64-bit shell and execute the following command
```
pacman -Syu
```
Run the following commands to install all pre-requisites
```
pacman -S mingw-w64-x86_64-gcc
pacman -S mingw-w64-x86_64-cmake
pacman -S make
pacman -S git
```
Clone the repo. Create a directory where the binaries will be built. The typical choice is a subfolder build. In this case, the commands are:
```
mkdir build
cd build
cmake ..
```
Build the program using command and run example
```
make
./bin/main
```
# How to compile and run on MacOS

Install the Mac terminal command line functions if needed (type git at the command line to trigger the install). Then install home-brew if not already present:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```
Install pre-requisites cmake library using Homebrew:
```
brew install cmake
```
Compile and run using dollowing commands (from the project directory)
```
mkdir build
cd build
cmake ..
make 
./bin/main
```
