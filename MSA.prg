define/parameter p1  CHAR 
define/parameter p2  N

set/form I1
!write/keyw i/I/1/1
def/loca i/I/1/1 0

!crea/disp
!load/itt neg
load/imag {p1}
clear/chan 2

!read/desc {p1} START >params.dat
!read/desc {p1} STEP >>params.dat


write/out "Find the galaxy center:"
get/curs tab_center.tbl
read/tabl tab_center.tbl :X_coordpix :Y_coordpix >tab_center.dat
$ rm tab_center.tbl

DO i = 1 {p2} 1
   write/out Sample several points on the {i} arm
   get/curs tab_arm{i}.tbl
   read/table tab_arm{i}.tbl :X_coordpix :Y_coordpix >arm{i}.dat
   $ rm tab_arm{i}.tbl
ENDDO

write/out "Computing..."
$ ./MSA.py {p1} {p2}
!bye