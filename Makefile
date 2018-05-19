1D_FEM.o : 1D_FEM.c  
	mpicc 1D_FEM.c -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -L /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib/  -Wl,-rpath,/home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib  -lpetsc -g

all: 1D_FEM.c 
	mpicc 1D_FEM.c -I /home/chen/Desktop/software/petsc/petsc-3.5.4/include -I /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/include -I /usr/lib/openmpi/include -L /home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib/  -Wl,-rpath,/home/chen/Desktop/software/petsc/petsc-3.5.4/x86_64/lib  -lpetsc -g



run:
	./a.out
	 

clean :
	rm *.o a.out
