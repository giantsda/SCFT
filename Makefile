nrutil.o: nrutil.c
	gcc nrutil.c -c -g

gaussj.o: gaussj.c  
	gcc gaussj.c  -c -g  

adm.o: adm.c    
	gcc adm.c -c -g


1D_FEM.o: 1D_FEM.c   
	gcc 1D_FEM.c -c  -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -g

all: 1D_FEM.o adm.o gaussj.o nrutil.o 
	mpicc 1D_FEM.o adm.o gaussj.o nrutil.o -o a.out -Wall -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -L /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib/  -Wl,-rpath,/home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib  -lpetsc -g -lm   -L /home/chen/Desktop/projects/shared_NRC -lhaha -Wl,-rpath,/home/chen/Desktop/projects/shared_NRC


run: 1D_FEM.o adm.o gaussj.o nrutil.o 
	mpicc 1D_FEM.o adm.o gaussj.o nrutil.o -o a.out -Wall -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -L /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib/  -Wl,-rpath,/home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib  -lpetsc -g -lm   -L /home/chen/Desktop/projects/shared_NRC -lhaha -Wl,-rpath,/home/chen/Desktop/projects/shared_NRC 
	./a.out
	 

clean :
	rm *.o a.out

 
