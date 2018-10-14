nrutil.o: nrutil.c
	gcc nrutil.c -c -g

gaussj.o: gaussj.c  
	gcc gaussj.c  -c -g  

adm.o: adm.c    
	gcc adm.c -c -g

broydn.o: broydn.c    
	gcc broydn.c -c -g

1D_FEM.o: 1D_FEM.c   
	gcc 1D_FEM.c -c  -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -g

1d_fem: 1D_FEM.o adm.o gaussj.o nrutil.o broydn.o  
	mpicc 1D_FEM.o adm.o gaussj.o nrutil.o broydn.o  -o 1d_fem -Wall -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -L /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib/  -Wl,-rpath,/home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib  -lpetsc -g -lm   -L . -lNR_C -Wl,-rpath,.

.PHONY: all

all:1d_fem


.PHONY: run

run:  
	make all
	./1d_fem
	 
.PHONY: clean

clean :
	rm *.o 1d_fem solution_store.txt

 
