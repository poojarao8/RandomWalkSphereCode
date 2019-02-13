##g++ -g test.cpp stl.cpp inv.cpp -o test.o -std=gnu++11 
##g++ main_parallel.cpp stl_parallel.cpp inv_parallel.cpp -o test.o -std=gnu++11

mpiCC main.cpp geometry.cpp randNumGen.cpp vecMatOps.cpp wlkrClass.cpp -o main.o -std=gnu++11


