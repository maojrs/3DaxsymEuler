! Module for aux averaged time dependant EOS parameters in aux arrays
MODULE auxmodule
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(:,:,:), Allocatable    :: ustar_array 
    CONTAINS
        
    SUBROUTINE init_auxmodule(mx,my,mbc,maux,aux)
        INTEGER :: mx,my,mbc,maux,i,j
        REAL(kind=8) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        Allocate(ustar_array(1-mbc:mx+mbc,1-mbc:my+mbc,1:2))
    END SUBROUTINE
    
END MODULE auxmodule