*     PROGRAM AAA9
********************
      use mpi
********************
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      PARAMETER (TOTAL= 500 *1000000)
********************
      character*10 M_filename
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
      integer::M_flag1,M_flag2
      INTEGER(KIND=MPI_OFFSET_KIND)::M_offset,M_shared_offset,M_offset_temp
********************
      COMMON /BIG/ F(TOTAL)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV,LNN
      COMMON /VALDIM/ FINDIM,TOT
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
**********************
      M_flag1 = 1
      M_flag2 = 1
      call mpi_init(M_ierr)
      call mpi_comm_rank(mpi_comm_world,M_myrank,M_ierr)
      call mpi_comm_size(mpi_comm_world,M_nproc,M_ierr)

********************
      IOUT=0
      ORDI=2
      MIL=1000000
      LIMV=MIL
      LNN=4*MIL
      TOT=TOTAL
      T1=TOT/1000000
      T2=0
      T3=0
********************
      if(M_myrank == 0) then
********************
      WRITE(IOUT,203)
203   FORMAT(10X,'SHELL MODEL CODE ANTOINE',/,10X,'by Etienne CAURIER',//,
     * 'THE USE OF THIS CODE IS AUTHORIZED UPON ACKNOWLEDGEMENT OF THE INTELLECTUAL PROPERTY ',/,
     * 'BY QUOTING THE FOLLOWING REFERENCES:',//,
     * 'E. CAURIER, shell model code ANTOINE,',/,
     * 'IRES, STRASBOURG 1989-2004',//,
     * 'E. CAURIER, F. NOWACKI',/,
     * 'Acta Physica Polonica 30 (1999) 705',//,
     * 'E. CAURIER, G. MARTINEZ-PINEDO, F. NOWACKI, A. POVES, A. P. ZUKER',/,
     * 'Reviews of Modern Physics 77, No 2, April 2005',//)
      WRITE(IOUT,202) T1,T2,T3
202   FORMAT(/,1X,'SIZE OF ARRAY F     =     ',I4,1X,I3.3,1X,I3.3,/)
      CALL ACAS (F)
********************
      else  
986         call mpi_recv(M_flag2,1,MPI_INTEGER,0,95,mpi_comm_world,M_status,M_ierr)
            if (M_flag2==2) then
                  go to 989
            end if
            call mpi_recv(M_filename,10,MPI_CHARACTER,0,97,mpi_comm_world,M_status,M_ierr)
            call mpi_file_open(mpi_comm_world,'mat/'//M_filename,MPI_MODE_RDONLY,MPI_INFO_NULL,M_myfile,M_ierr)
*            write(*,*) 'recv98'
987         call mpi_recv(M_flag1,1,MPI_INTEGER,0,96,mpi_comm_world,M_status,M_ierr)
            if (M_flag1==2) then
                  go to 988
            end if
            call mpi_recv(Mi_buf,2,MPI_INTEGER,0,98,mpi_comm_world,M_status,M_ierr)
            M_ZSIZE=Mi_buf(1)
            M_shared_offset=Mi_buf(2)
            call mpi_file_seek(M_myfile,((M_myrank-1)*(M_ZSIZE/(M_nproc-1))*2*sizeof(mpi_double_precision))+M_shared_offset,MPI_SEEK_SET,M_status,M_ierr)
            call mpi_file_read(M_myfile,Md_buf,M_ZSIZE/(M_nproc-1),MPI_DOUBLE_PRECISION,M_status,M_ierr)
            do i = 1,M_ZSIZE/(M_nproc-1)
*            write(*,*) 'send99:',Md_buf(i),M_myrank
            end do
            call mpi_send(Md_buf,M_ZSIZE/(M_nproc-1),MPI_DOUBLE_PRECISION,0,99,mpi_comm_world,M_status,M_ierr)
            GO TO 987
988         call mpi_file_close(M_myfile,M_ierr)
                  go to 986
********************
********************
989         write (*,*) '********************MPIFINISHNON0'
            call mpi_finalize(M_ierr)
            
      end if
********************
      STOP
      END
