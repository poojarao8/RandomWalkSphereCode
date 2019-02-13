##g++ -g test.cpp stl.cpp inv.cpp -o test.o -std=gnu++11 
g++ -pg main_serial.cpp stl_serial.cpp inv_serial.cpp -o main_serial.o -std=gnu++11

./main_serial.o
