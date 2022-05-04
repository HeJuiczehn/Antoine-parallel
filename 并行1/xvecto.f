      SUBROUTINE OVERL (V1,V2,YOVER,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V1(*),V2(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      YOVER=0.0D0
      DO 2  N=1,NDIM
2     YOVER=YOVER+V1(N)*V2(N)
      ELSE
      CALL SOVERL (V1,V2,YOVER,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SOVERL (U1,U2,YOVER,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U1(*),U2(*)
      COMMON IOUT
      YOVER=0.0D0
      DO 1  N=1,NDIM
1     YOVER=YOVER+DBLE(U1(N))*DBLE(U2(N))
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE OVERM (V1,V2,YOVER,KMIN,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V1(*),V2(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      YOVER=0.0D0
      DO 2  N=KMIN,NDIM
2     YOVER=YOVER+V1(N)*V2(N)
      ELSE
      CALL SOVERM (V1,V2,YOVER,KMIN,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SOVERM (U1,U2,YOVER,KMIN,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U1(*),U2(*)
      COMMON IOUT
      YOVER=0.0D0
      DO 1  N=KMIN,NDIM
1     YOVER=YOVER+DBLE(U1(N))*DBLE(U2(N))
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VPSCAL (V,YSCAL,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      YSCAL=0.0D0
      DO 2  N=1,NDIM
2     YSCAL=YSCAL+V(N)*V(N)
      ELSE
      CALL SPSCAL (V,YSCAL,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SPSCAL (U,YSCAL,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      YSCAL=0.0D0
      DO 1  N=1,NDIM
1     YSCAL=YSCAL+DBLE(U(N))*DBLE(U(N))
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VMULT (V,ZCOEF,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 2  N=1,NDIM
2     V(N)=V(N)*ZCOEF
      ELSE
      CALL SMULT (V,ZCOEF,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SMULT (U,ZCOEF,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*)
      COMMON IOUT
      UCOEF=REAL(ZCOEF)
      DO 1  N=1,NDIM
1     U(N)=U(N)*UCOEF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VLINE (V1,V2,ZCOEF,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V1(*),V2(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 2  N=1,NDIM
2     V1(N)=V1(N)+V2(N)*ZCOEF
      ELSE
      CALL SLINE (V1,V2,ZCOEF,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SLINE (U1,U2,ZCOEF,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U1(*),U2(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      UCOEF=REAL(ZCOEF)
      DO 1  N=1,NDIM
1     U1(N)=U1(N)+U2(N)*UCOEF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VADD (V1,V2,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V1(*),V2(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 2  N=1,NDIM
2     V1(N)=V1(N)+V2(N)
      ELSE
      CALL SADD (V1,V2,NDIM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SADD (U1,U2,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U1(*),U2(*)
      COMMON IOUT
      DO 1  N=1,NDIM
1     U1(N)=U1(N)+U2(N)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE ORTHOG (F,FF,YOVER,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*),FF(*)
      COMMON IOUT
      CALL OVERL (F,FF,YOVER,NDIM,KPREC)
      CALL VLINE (F,FF,-YOVER,NDIM,KPREC)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VNORM (F,ZN,NDIM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      CALL VPSCAL (F,ZN,NDIM,KPREC)
      IF(ZN.GT.1.0D-10) THEN
      ZC=1.0D0/DSQRT(ZN)
      CALL VMULT (F,ZC,NDIM,KPREC)
      ELSE
      WRITE(IOUT,*) 'NORM OF THE VECTOR=0'
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE INIFIX (A,N,Q)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION A(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      DO 1 I=1,N
1     A(I)=Q
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE IPRIN (A,N)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION A(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      WRITE(IOUT,*) (A(K),K=1,N)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIMPLE (U,V,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*),V(*)
      COMMON IOUT
      DO 1 N=1,NDIM
1     U(N)=REAL(V(N))
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DOUBLE (V,U,NDIM)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*),V(*)
      COMMON IOUT
      DO 1 N=1,NDIM
1     V(N)=DBLE(U(N))
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VINI0 (V,N,KPREC)
      DOUBLE PRECISION  V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 1 I=1,N
1     V(I)=0.0D0
      ELSE
      CALL SINI0 (V,N)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SINI0 (U,N)
      DIMENSION U(*)
      COMMON IOUT
      DO 1 I=1,N
1     U(I)=0.0
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VINIZ (V,N,Z,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 1 I=1,N
1     V(I)=Z
      ELSE
      CALL SINIZ (V,N)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SINIZ (U,N,Z)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*)
      COMMON IOUT
      UU=REAL (Z)
      DO 1 I=1,N
1     U(I)=UU
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VPRIN (V,N,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      WRITE(IOUT,200) (V(K),K=1,N)
200   FORMAT(6F11.6)
      ELSE
      CALL SPRIN (V,N)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SPRIN (U,N)
      IMPLICIT INTEGER (A-T)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION U(*)
      COMMON IOUT
      WRITE(IOUT,200) (U(K),K=1,N)
200   FORMAT(6F11.6)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VMOVE (W,V,N,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION W(*),V(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      DO 1 I=1,N
1     W(I)=V(I)
      ELSE
      CALL SMOVE (W,V,N)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SMOVE (W,V,N)
      IMPLICIT INTEGER (A-U)
      DIMENSION W(*),V(*)
      COMMON IOUT
      DO 1 I=1,N
1     W(I)=V(I)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE IMOVE (A,B,N)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION A(*),B(*)
      COMMON IOUT
      DO 1 I=1,N
1     A(I)=B(I)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE IMAXI (MAXI,LIST,K)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      COMMON IOUT
      DIMENSION LIST(*)
      MAXI=LIST(1)
      DO 1 I=2,K
      MAXI=MAX(MAXI,LIST(I))
1     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE VMAXI (MAXI,Z,K,IPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT
      IF(IPREC.EQ.2) THEN
      ZMAX=Z(1)
      MAXI=1
      DO 1 I=2,K
      IF(Z(I).LT.ZMAX) GO TO 1
      ZMAX=Z(I)
      MAXI=I
1     CONTINUE
      ELSE
      CALL SMAXI (MAXI,Z,K)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SMAXI (MAXI,Z,K)
      IMPLICIT INTEGER (A-U)
***   IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT
      ZMAX=Z(1)
      MAXI=1
      DO 1 I=2,K
      IF(Z(I).LT.ZMAX) GO TO 1
      ZMAX=Z(I)
      MAXI=I
1     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE LECFIX (FLLEC,F,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      READ (FLLEC,END=99) (F(K),K=1,MIN(LIMV,NTERM))
      IF(NTERM.GT.LIMV) THEN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      READ (FLLEC,END=99) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
99    WRITE(IOUT,*) 'END READING IN LECFIX  ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECFIX_2 (FLLEC,F,NTERM)
      use mpi
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
********************
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
********************
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
**********************
****************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Mi_buf,MIN(LIMV,NTERM),MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      do K=1,MIN(LIMV,NTERM)
            F(K)=Mi_buf(K)
      end do
****************************
*      READ (FLLEC,END=99) (F(K),K=1,MIN(LIMV,NTERM))
      IF(NTERM.GT.LIMV) THEN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
****************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Mi_buf,NMAX-NMIN+1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      do K=NMIN,NMAX
            F(K)=Mi_buf(K-NMIN+1)
      end do
****************************
*      READ (FLLEC,END=99) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
99    WRITE(IOUT,*) 'END READING IN LECFIX  ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECVEC (FLLEC,Z,NTERM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      CALL LECDBL (FLLEC,Z,NTERM)
      ELSE
      CALL LECSMP (FLLEC,Z,NTERM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE LECVEC_2 (FLLEC,Z,NTERM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      CALL LECDBL_2 (FLLEC,Z,NTERM)
      ELSE
      CALL LECSMP_2 (FLLEC,Z,NTERM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE LECSMP (FLLEC,F,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      REAL*4    F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      READ (FLLEC,END=99) (F(K),K=1,MIN(LIMV,NTERM))
      IF(NTERM.GT.LIMV) THEN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      READ (FLLEC,END=99) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
99    WRITE(IOUT,*) 'END OF READING IN LECSMP      ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECSMP_2 (FLLEC,F,NTERM)
      use mpi
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      REAL*4    F(*)
********************
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
********************
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
**********************
****************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Mr_buf,MIN(LIMV,NTERM),MPI_REAL4,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      do K=1,MIN(LIMV,NTERM)
            F(K)=Mr_buf(K)
      end do
****************************
*      READ (FLLEC,END=99) (F(K),K=1,MIN(LIMV,NTERM))
      IF(NTERM.GT.LIMV) THEN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
****************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Mr_buf,NMAX-NMIN+1,MPI_REAL4,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      do K=NMIN,NMAX
            F(K)=Mr_buf(K-NMIN+1)
      end do
****************************
*      READ (FLLEC,END=99) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
99    WRITE(IOUT,*) 'END OF READING IN LECSMP      ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECDBL (FLLEC,Z,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      READ (FLLEC,END=99) (Z(K),K=1,MIN(LIMV,NTERM))
2     FORMAT(7F10.5)
      IF(NTERM.LE.LIMV) RETURN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      READ (FLLEC,END=99) (Z(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      RETURN
99    WRITE(IOUT,*) 'END OF READING IN LECDBL  ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECDBL_2 (FLLEC,Z,NTERM)
      use mpi
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
********************
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
      integer*4::M_s
      INTEGER (KIND=MPI_OFFSET_KIND)::M_offset,M_shared_offset
********************
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
**********************
***************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)

      call mpi_file_get_position(M_myfile,M_shared_offset,M_ierr)
      if (M_nproc > 1) then
            Mi_buf(1)=MIN(LIMV,NTERM)
            Mi_buf(2)=M_shared_offset
            do index=1,M_nproc-1
                  call mpi_send(Mi_buf,2,MPI_INTEGER,index,98,mpi_comm_world,M_status,M_ierr)
            end do
            do index=1,M_nproc-1
                  call mpi_recv(Md_buf,MIN(LIMV,NTERM)/(M_nproc-1),MPI_DOUBLE_PRECISION,index,99,mpi_comm_world,M_status,M_ierr)
                  do index2=1,MIN(LIMV,NTERM)/(M_nproc-1)
                        Z((index-1)*(MIN(LIMV,NTERM)/(M_nproc-1))+index2)=Md_buf(index2)
*                        write(*,*) Md_buf(index2)
                  end do
            end do
      end if
      
      if (M_nproc > 1) then
            M_offset = (M_nproc-1)*(MIN(LIMV,NTERM)/(M_nproc-1))*2*sizeof(mpi_double_precision)+M_shared_offset
            call MPI_FILE_SEEK(M_myfile,M_offset,MPI_SEEK_SET,M_ierr)
            if (mod(MIN(LIMV,NTERM),(M_nproc-1)) > 0) then

                  call MPI_FILE_READ(M_myfile,Md_buf,mod(MIN(LIMV,NTERM),(M_nproc-1)),MPI_DOUBLE_PRECISION,M_status,M_ierr)

                  do index2=1,mod(MIN(LIMV,NTERM),(M_nproc-1))
                        Z((M_nproc-1)*(MIN(LIMV,NTERM)/(M_nproc-1))+index2)=Md_buf(index2)
*                        write(*,*) 'redundancy:',Md_buf(index2)
                  end do
            end if
            
      else
            call MPI_FILE_READ(M_myfile,Md_buf,MIN(LIMV,NTERM),MPI_DOUBLE_PRECISION,M_status,M_ierr)

            do index2=1,MIN(LIMV,NTERM)
                  Z(index2)=Md_buf(index2)
            end do
      end if
                  
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)

      
***************************
*      READ (FLLEC,END=99) (Z(K),K=1,MIN(LIMV,NTERM))
2     FORMAT(7F10.5)
      IF(NTERM.LE.LIMV) RETURN
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
****************************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Md_buf,NMAX-NMIN+1,MPI_DOUBLE_PRECISION,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      do K=NMIN,NMAX
            Z(K)=Md_buf(K-NMIN+1)
*            write(*,*)Z(K)
      end do
****************************
*      READ (FLLEC,END=99) (Z(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      RETURN
99    WRITE(IOUT,*) 'END OF READING IN LECDBL  ',NTERM
      STOP
      END
*
************************************************************************
*
      SUBROUTINE WRIFIX (FLWRI,F,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      IF(NTERM.LE.LIMV) THEN
      WRITE (FLWRI) (F(K),K=1,NTERM)
      ELSE
      WRITE (FLWRI) (F(K),K=1,LIMV)
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      WRITE (FLWRI) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
      END

*
************************************************************************
*
      SUBROUTINE WRIVEC (FLWRI,Z,NTERM,KPREC)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT
      IF(KPREC.EQ.2) THEN
      CALL WRIDBL (FLWRI,Z,NTERM)
      ELSE
      CALL WRISMP (FLWRI,Z,NTERM)
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE WRISMP (FLWRI,F,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      REAL*4 F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      IF(NTERM.LE.LIMV) THEN
      WRITE (FLWRI) (F(K),K=1,NTERM)
      ELSE
      WRITE (FLWRI) (F(K),K=1,LIMV)
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      WRITE (FLWRI) (F(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE WRIDBL (FLWRI,Z,NTERM)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION Z(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /LIMIT/ LIMV
      IF(NTERM.LE.LIMV) THEN
      WRITE (FLWRI) (Z(K),K=1,NTERM)
      ELSE
      WRITE (FLWRI) (Z(K),K=1,LIMV)
      NMAX=LIMV
1     NMIN=NMAX+1
      NMAX=MIN(NTERM,NMAX+LIMV)
      WRITE (FLWRI) (Z(K),K=NMIN,NMAX)
      IF(NMAX.LT.NTERM) GO TO 1
      END IF
      RETURN
      END
