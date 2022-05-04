      SUBROUTINE PPMAT (VFIN,VINI,F)
********************
      use mpi
********************
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
********************
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
********************
      REAL*4  US,USEC
      CHARACTER*10 NOM,FIL
      DIMENSION F(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB,DO,O
      COMMON /BLOC/ NBL
      COMMON /BLOCA/ TI(585,4),TA(585,4)
      COMMON /FONCT/ JUMP,TFIX,NDIM
      COMMON /KAS/ NF,HJT
      COMMON /VALDIM/ FINDIM,TOT,FINFON,FIN0
      COMMON /VALSL/ NSL
      COMMON /TICAL/ USEC(9)
      COMMON /NAME/ NOM,FIL
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
*      write(*,*) "PPMAT BEGIN"
**********************
      US=SECNDS (0.0)
      FIL='ZRF       '
      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
     *       STATUS='OLD')
      RF=FINDIM+1
      RMIN=RF+NSL
      RMAX=RMIN+NSL
      FIL=NOM
      LL=3*NSL
      FIN0=FINDIM+LL
      CALL LECFIX (1,F(RF),LL)
      CLOSE (1,STATUS='KEEP')
      NBI=FIN0+1
      IDIV=0
      ISTO=0
      WRITE(FIL(3:3),101) 1
101   FORMAT(I1)
      WRITE(FIL(4:5),100) IDIV
100   FORMAT(I2)
      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
     *       STATUS='OLD')
1     CONTINUE
      ISTO=ISTO+1
      IF(ISTO.GT.30) THEN
      CLOSE (1,STATUS='KEEP')
      ISTO=1
      IDIV=IDIV+1
      WRITE(FIL(4:5),100) IDIV
      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
     *       STATUS='OLD')
      END IF
      CALL LECNN2 (IFINAL,NTERM,F)
      IF(NTERM.EQ.0) GO TO 2
      IF(IPREC.EQ.2) THEN
      CALL DIANN2 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),VFIN,VINI)
      ELSE
      CALL SIANN2 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),VFIN,VINI)
      END IF
      IF(IFIN.NE.IFINAL) GO TO 1
2     CLOSE (1,STATUS='KEEP')
      US=SECNDS (US)
      USEC(1)=USEC(1)+US
**    WRITE(IOUT,*) 'FIN DIANN2'
***********************************
      US=SECNDS (0.0)
      FIL=NOM
      WRITE(FIL(3:3),101) 2
      DO 3 NB=1,NBL
      WRITE(FIL(4:5),100) NB
*****************
      if (M_nproc > 1) then
            do index=1,M_nproc-1
                  call mpi_send(1,1,MPI_INTEGER,index,95,mpi_comm_world,M_status,M_ierr)
            end do
      end if
      if (M_nproc > 1) then
            do index=1,M_nproc-1
                  call mpi_send(FIL,10,MPI_CHARACTER,index,97,mpi_comm_world,M_status,M_ierr)
            end do
      end if
      call mpi_file_open(mpi_comm_world,'mat/'//FIL,MPI_MODE_RDONLY,MPI_INFO_NULL,M_myfile,M_ierr)

*****************
*      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
*     *       STATUS='OLD')
      DO 4 Q2=TI(NB,2),TA(NB,2)
      CALL BORQ (2,1,NB,Q2,Q1MIN,Q1MAX)
*****************
*      write(*,*) "lecnn2 begin"
5     CALL LECNN2_2 (IFINAL,NTERM,F)
*      write(*,*) "lecnn2 end"
*****************
      IF(NTERM.EQ.0) GO TO 4
*      write(*,*) "+++++++++++++mark1"
      IF(IPREC.EQ.2) THEN
*      write(*,*) "+++++++++++++mark2"
*      write(*,*) Q1MIN,Q1MAX,F(NBI),F(RF),F(DO),F(O)
      CALL DIAPP2 (Q1MIN,Q1MAX,F(NBI),F(RF),F(DO),F(O),VFIN,VINI)
*      write(*,*) "+++++++++++++mark3"
      ELSE
      CALL SIAPP2 (Q1MIN,Q1MAX,F(NBI),F(RF),F(DO),F(O),VFIN,VINI)
*      write(*,*) "+++++++++++++mark4"
      END IF
      IF(IFIN.NE.IFINAL) GO TO 5
*      write(*,*) "+++++++++++++mark5"
4     CONTINUE
*****************
      if (M_nproc > 1) then
            do index=1,M_nproc-1
                  call mpi_send(2,1,MPI_INTEGER,index,96,mpi_comm_world,M_status,M_ierr)
            end do
      end if
*      write (*,*) 'file closing'
      call mpi_file_close(M_myfile,M_ierr)
*****************
*     CLOSE (1,STATUS='KEEP')
3     CONTINUE
      US=SECNDS (US)
      USEC(2)=USEC(2)+US
**    WRITE(IOUT,*) 'FIN DIAPP2'
*****************************************
      IF(HJT.NE.0) RETURN
      DO 10 PN=1,2
      US=SECNDS (0.0)
      FIL=NOM
      WRITE(FIL(3:3),101) PN+2
      IF(PN.EQ.1) THEN
      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
     *       STATUS='UNKNOWN')
      END IF
      DO 6 NB=1,NBL
      IF(PN.EQ.2) THEN
      WRITE(FIL(4:5),100) NB
      OPEN (UNIT=1,FILE='mat/'//FIL,FORM='UNFORMATTED',
     *       STATUS='UNKNOWN')
      END IF
      READ (1,END=99) IDEB,IFIN,NTERM
***   WRITE(IOUT,*) 'NB',NB,IDEB,IFIN,NTERM
      LL=IFIN-IDEB+1
      READ (1,END=99) (F(FIN0+K),K=1,LL)
      IF(NTERM.EQ.0) GO TO 6
      FIN1=FIN0+LL
      FIN1=FIN1+MOD(FIN1,ORDI)
      DO=FIN1+1
      O=DO+NTERM*IPREC
      FIN1=FIN1+(IPREC+2)*NTERM
      CALL LECVEC (1,F(DO),NTERM,IPREC)
      CALL LECFIX (1,F(O),NTERM*2)
      ONE=FIN1+1
      READ (1,END=99) NONE,KV,(F(FIN1+K),K=1,KV)
      FIN1=FIN1+KV
      FIN1=FIN1+MOD(FIN1,ORDI)
      DONE=FIN1+1
      CALL TTDIM (FIN1+NONE*IPREC,TOT,'PPMATD')
      CALL LECVEC (1,F(DONE),NONE,IPREC)
      IF(IPREC.EQ.2) THEN
      IF(PN.EQ.1) THEN
      IF(TFIX.EQ.0) THEN
      CALL DIANN1 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),F(DONE),F(ONE),
     *             VFIN,VINI)
      ELSE
      CALL DIATT1 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),F(DONE),F(ONE),
     *             VFIN,VINI)
      END IF
      ELSE
      CALL DIAPP1 (F(NBI),F(RF),F(DO),F(O),F(DONE),F(ONE),VFIN,VINI)
      CLOSE (1,STATUS='KEEP')
      END IF
      ELSE
      IF(PN.EQ.1) THEN
      IF(TFIX.EQ.0) THEN
      CALL SIANN1 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),F(DONE),F(ONE),
     *             VFIN,VINI)
      ELSE
      CALL SIATT1 (F(NBI),F(RMIN),F(RMAX),F(DO),F(O),F(DONE),F(ONE),
     *             VFIN,VINI)
      END IF
      ELSE
      CALL SIAPP1 (F(NBI),F(RF),F(DO),F(O),F(DONE),F(ONE),VFIN,VINI)
      CLOSE (1,STATUS='KEEP')
      END IF
      END IF
6     CONTINUE
      IF(PN.EQ.1) CLOSE (1,STATUS='KEEP')
      US=SECNDS (US)
      USEC(PN+2)=USEC(PN+2)+US
**    WRITE(IOUT,*) 'FIN DIAPP1',PN
10    CONTINUE
      RETURN
99    WRITE(IOUT,*) 'END OF READING IN PPMAT  FILE  ',FIL
      STOP
      END
*
************************************************************************
*
      SUBROUTINE LECNN2 (IFINAL,NTERM,F)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB,DO,O
      COMMON /VALDIM/ FINDIM,TOT,FINFON,FIN0
      READ  (1) IDEB,IFIN,IFINAL,NTERM
      IF(NTERM.EQ.0) RETURN
      LL=IFIN-IDEB+1
      READ (1) (F(FIN0+K),K=1,LL)
      FIN1=FIN0+LL
      FIN1=FIN1+MOD(FIN1,ORDI)
      DO=FIN1+1
      O=DO+NTERM*IPREC
      CALL TTDIM (FIN1+(IPREC+1)*NTERM,TOT,'PPMATA')
      CALL LECVEC (1,F(DO),NTERM,IPREC)
      CALL LECFIX (1,F(O),NTERM)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE LECNN2_2 (IFINAL,NTERM,F)
********************
      use mpi
********************
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
********************
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
********************
      DIMENSION F(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB,DO,O
      COMMON /VALDIM/ FINDIM,TOT,FINFON,FIN0
********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPISTATUS/ M_status
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc
      COMMON /MPIFILE/ M_myfile
********************
*******************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
* unformatted file has one redundant int data which record the length
      call MPI_FILE_READ(M_myfile,Mi_buf,4,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      IDEB=Mi_buf(1)
      IFIN=Mi_buf(2)
      IFINAL=Mi_buf(3)
      NTERM=Mi_buf(4)
*      write(*,*) IDEB,IFIN,IFINAL,NTERM
*      write(*,*) F(DO),F(O)

*******************
*      READ  (1) IDEB,IFIN,IFINAL,NTERM
*      write(*,*) IDEB,IFIN,IFINAL,NTERM
      IF(NTERM.EQ.0) RETURN
      if (M_nproc > 1) then
            do index=1,M_nproc-1
                  call mpi_send(1,1,MPI_INTEGER,index,96,mpi_comm_world,M_status,M_ierr)
            end do
      end if
      LL=IFIN-IDEB+1
*******************
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,Mi_buf,LL,MPI_INTEGER,M_status,M_ierr)
      call MPI_FILE_READ(M_myfile,M_s,1,MPI_INTEGER,M_status,M_ierr)
*      write(*,*) LL
      do K=1,LL
            F(FIN0+K)=Mi_buf(K)
*            write(*,*) Mi_buf(K)
      end do

*******************
*      READ (1) (F(FIN0+K),K=1,LL)
      FIN1=FIN0+LL
      FIN1=FIN1+MOD(FIN1,ORDI)
      DO=FIN1+1
      O=DO+NTERM*IPREC
      CALL TTDIM (FIN1+(IPREC+1)*NTERM,TOT,'PPMATA')
*      write(*,*) "F(DO)andF(O)before lecvec and lecfix",F(DO),F(O)
      CALL LECVEC_2 (1,F(DO),NTERM,IPREC)
*      write(*,*) "F(DO)andF(O)after lecvec",F(DO),F(O)
      CALL LECFIX_2 (1,F(O),NTERM)
*      write(*,*) "F(DO)andF(O)after lecfix",F(DO),F(O)
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DIANN2 (NBI,RMIN,RMAX,ZO,O,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION ZO(*),O(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN
      K=0
      DO 2 I1=IDEB,IFIN
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 2 Q=1,NBI(I1)
      K=K+1
      DEL=O(K)
      ZZ=ZO(K)
      DO 2 N=NMIN,NMAX
      NN=N+DEL
      VFIN(N)=VFIN(N)+VINI(NN)*ZZ
      VFIN(NN)=VFIN(NN)+VINI(N)*ZZ
***   WRITE(IOUT,*) 'NN2',N,NN,ZZ
2     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DIAPP2 (Q1MIN,Q1MAX,NBI,RF,ZO,O,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(*),ZO(*),NBI(IDEB:*),RF(*)
      DIMENSION VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      DO 10 Q1=Q1MIN,Q1MAX
      DO 1 I1=DEB(Q1,1),FIN(Q1,1)
      RI1=RF(I1)
      K=0
      DO 1 I2=IDEB,IFIN
      N=RI1+I2
      DO 1 Q=1,NBI(I2)
      K=K+1
      NN=N+O(K)
      ZZ=ZO(K)
      VFIN(N)=VFIN(N)+VINI(NN)*ZZ
      VFIN(NN)=VFIN(NN)+VINI(N)*ZZ
**    WRITE(IOUT,*) 'PP2',N,NN,ZZ
1     CONTINUE
10    CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DIANN1 (NBI,RMIN,RMAX,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN
      K=0
      DO 4 I1=IDEB,IFIN
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 4 Q=1,NBI(I1)
      K=K+1
      DEL=O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=ONE(T)
      DO 5 N=NMIN,NMAX
      NN=N+DEL
      R=R+1
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
5     CONTINUE
4     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DIATT1 (NBI,RMIN,RMAX,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      COMMON /BLOCA/ TI(585,4),TA(585,4)
      K=0
      DO 4 Q1=TI(NB,1),TA(NB,1)
      CALL BORQ (1,1,NB,Q1,Q2MIN,Q2MAX)
      RR=0
      IF(Q2MIN.NE.TI(NB,2)) THEN
      DO 2 Q2=TI(NB,2),Q2MIN-1
2     RR=RR+SIZ(Q2,2)
      END IF
      DO 4 I1=DEB(Q1,1),FIN(Q1,1)
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 4 Q=1,NBI(I1)
      K=K+1
      DEL=O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=RR+ONE(T)
      DO 5 N=NMIN,NMAX
      NN=N+DEL
      R=R+1
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
5     CONTINUE
4     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE DIAPP1 (NBI,RF,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RF(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      COMMON /BLOCA/ TI(585,4),TA(585,4)
      COMMON /FONCT/ JUMP,TFIX,NDIM
      RR=0
      DO 1 Q1=TI(NB,1),TA(NB,1)
      KINI=0
      CALL BORQ (1,1,NB,Q1,Q2MIN,Q2MAX)
      IF(TFIX.NE.0) THEN
      DO 7 I2=DEB(TI(NB,2),2),FIN(Q2MIN-1,2)
      KINI=KINI+NBI(I2)
7     CONTINUE
      END IF
      DO 8 I1=DEB(Q1,1),FIN(Q1,1)
      RR=RR+1
      RI1=RF(I1)
      K=KINI
      DO 8 I2=DEB(Q2MIN,2),FIN(Q2MAX,2)
      N=RI1+I2
      DO 8 Q=1,NBI(I2)
      K=K+1
      NN=N+O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=ONE(T)+RR
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
8     CONTINUE
1     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIANN2 (NBI,RMIN,RMAX,ZO,O,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
**    IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION ZO(*),O(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN
      K=0
      DO 2 I1=IDEB,IFIN
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 2 Q=1,NBI(I1)
      K=K+1
      DEL=O(K)
      ZZ=ZO(K)
      DO 2 N=NMIN,NMAX
      NN=N+DEL
      VFIN(N)=VFIN(N)+VINI(NN)*ZZ
      VFIN(NN)=VFIN(NN)+VINI(N)*ZZ
2     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIAPP2 (Q1MIN,Q1MAX,NBI,RF,ZO,O,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
**    IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(*),ZO(*),NBI(IDEB:*),RF(*)
      DIMENSION VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      DO 10 Q1=Q1MIN,Q1MAX
      DO 1 I1=DEB(Q1,1),FIN(Q1,1)
      RI1=RF(I1)
      K=0
      DO 1 I2=IDEB,IFIN
      N=RI1+I2
      DO 1 Q=1,NBI(I2)
      K=K+1
      NN=N+O(K)
      ZZ=ZO(K)
      VFIN(N)=VFIN(N)+VINI(NN)*ZZ
      VFIN(NN)=VFIN(NN)+VINI(N)*ZZ
1     CONTINUE
10    CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIANN1 (NBI,RMIN,RMAX,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
**    IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN
      K=0
      DO 4 I1=IDEB,IFIN
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 4 Q=1,NBI(I1)
      K=K+1
      DEL=O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=ONE(T)
      DO 5 N=NMIN,NMAX
      NN=N+DEL
      R=R+1
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
5     CONTINUE
4     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIATT1 (NBI,RMIN,RMAX,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
**    IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RMIN(*),RMAX(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      COMMON /BLOCA/ TI(585,4),TA(585,4)
      K=0
      DO 4 Q1=TI(NB,1),TA(NB,1)
      CALL BORQ (1,1,NB,Q1,Q2MIN,Q2MAX)
      RR=0
      IF(Q2MIN.NE.TI(NB,2)) THEN
      DO 2 Q2=TI(NB,2),Q2MIN-1
2     RR=RR+SIZ(Q2,2)
      END IF
      DO 4 I1=DEB(Q1,1),FIN(Q1,1)
      NMIN=RMIN(I1)
      NMAX=RMAX(I1)
      DO 4 Q=1,NBI(I1)
      K=K+1
      DEL=O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=RR+ONE(T)
      DO 5 N=NMIN,NMAX
      NN=N+DEL
      R=R+1
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
5     CONTINUE
4     CONTINUE
      RETURN
      END
*
************************************************************************
*
      SUBROUTINE SIAPP1 (NBI,RF,ZO,O,ZONE,ONE,VFIN,VINI)
      IMPLICIT INTEGER (A-U)
**    IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION O(2,*),ZO(*),NBI(IDEB:*),RF(*)
      DIMENSION ONE(*),ZONE(*),VINI(*),VFIN(*)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /AAA1/ IDEB,IFIN,NB
      COMMON /BORN/ DEB(888,4),FIN(888,4),SIZ(888,4)
      COMMON /BLOCA/ TI(585,4),TA(585,4)
      COMMON /FONCT/ JUMP,TFIX
      RR=0
      DO 1 Q1=TI(NB,1),TA(NB,1)
      KINI=0
      CALL BORQ (1,1,NB,Q1,Q2MIN,Q2MAX)
      IF(TFIX.NE.0) THEN
      DO 7 I2=DEB(TI(NB,2),2),FIN(Q2MIN-1,2)
      KINI=KINI+NBI(I2)
7     CONTINUE
      END IF
      DO 8 I1=DEB(Q1,1),FIN(Q1,1)
      RR=RR+1
      RI1=RF(I1)
      K=KINI
      DO 8 I2=DEB(Q2MIN,2),FIN(Q2MAX,2)
      N=RI1+I2
      DO 8 Q=1,NBI(I2)
      K=K+1
      NN=N+O(1,K)
      ZZ=ZO(K)
      T=O(2,K)
      R=ONE(T)+RR
      YY=ZZ+ZONE(R)
      VFIN(N)=VFIN(N)+VINI(NN)*YY
      VFIN(NN)=VFIN(NN)+VINI(N)*YY
8     CONTINUE
1     CONTINUE
      RETURN
      END
*
************************************************************************
*
