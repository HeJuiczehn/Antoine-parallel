      SUBROUTINE ACAS (F)
      use mpi
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      DIMENSION F(*)
      REAL*4 W
      CHARACTER*40 TEXT
      CHARACTER*20 NEWT
********************
      character*10 M_filename
      integer :: M_status(mpi_status_size)
      double precision,dimension(1024*1024) :: Md_buf
      real*4,dimension(1024*1024) :: Mr_buf
      integer,dimension(1024*1024) :: Mi_buf
      integer::M_flag1,M_flag2
      INTEGER (KIND=MPI_OFFSET_KIND)::M_offset,M_shared_offset
********************
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /BARAT/ ITEXT,TEXT(2),NEWT
      COMMON /FONCT/ A(13)
      COMMON /VALDIM/ FINDIM
**********************
      COMMON /MPIBUF/ Md_buf,Mr_buf,Mi_buf
      COMMON /MPIS/ M_ierr,M_myrank,M_nproc 
      COMMON /MPIFILE/ M_myfile
      COMMON /MPISTATUS/ M_status
**********************

      CALL ZCONST
1     READ (5,*,END=99) CAS,ITEXT,IPRI
      IF(CAS.EQ.0) GO TO 99
      IPREC=2
      IF(CAS.GT.100) THEN
      CAS=MOD(CAS,100)
      IPREC=1
      END IF
      A(13)=IPREC
      IF(ITEXT.NE.0) THEN
      READ (5,100,END=99) NEWT
100   FORMAT(A20)
      ELSE
      NEWT='                    '
      END IF
      WRITE(IOUT,200) CAS
200   FORMAT(/,2X,'*********** CAS=',I3,2X,'********************',/)
      IF(IPREC.EQ.1) WRITE(IOUT,*) 'CALCULATION WITH SIMPLE PRECISION'
      FINDIM=0
      W=SECNDS (0.0)
      IF(CAS.LE.40) THEN
      CALL ACHOI (F)
      ELSE
      CALL UTIL  (F)
      END IF
      W=SECNDS (W)
      WRITE(IOUT,201) W
201   FORMAT(/,2X,'TIME FOR THIS OPTION=',F15.3)
      WRITE(IOUT,202)
202   FORMAT(2X,20('***'))
      GO TO 1
99    WRITE(IOUT,*) 'NORMAL END OF THE JOB     '
**********
      write (*,*) '********************MPIFINISH0'
      if (M_nproc > 1) then
            do index=1,M_nproc-1
                  call mpi_send(2,1,MPI_INTEGER,index,95,mpi_comm_world,M_status,M_ierr)
            end do
      end if
      call mpi_finalize(M_ierr)
**********
      STOP
      END
*
************************************************************************
*
      SUBROUTINE ZCONST
      IMPLICIT INTEGER (A-U)
      IMPLICIT DOUBLE PRECISION (V-Z)
      COMMON IOUT,IPRI,ORDI,IPREC,CAS
      COMMON /CONST/ ZAC,ZACM,ZPI,ZEPS
      COMMON /CLOG/ ZLOG(100),ZGAMM2(100)
      COMMON /CHAT/ ZHAT(0:200)
      ZAC=SQRT(2.0D0)
      ZACM=1.0D0/ZAC
      ZPI=4.0D0*DATAN(1.0D0)
***   ZEPS=X02AJF(ZZ)
      ZEPS=1.0D-12
      ZLOG(1)=0.0D0
      ZLOG(2)=0.0D0
      WN=1.0D0
      DO 1 I=3,100
      WN=WN+1.0D0
1     ZLOG(I)=ZLOG(I-1)+DLOG(WN)
      zgamm2(1)=dsqrt(zpi)
      zgamm2(2)=1.0d0
      do jj=3,100
         zgamm2(jj)=(jj-2)*zgamm2(jj-2)/2.0d0
      end do
      DO 3 S=0,200
3     ZHAT(S)=DSQRT(DBLE(S+1))
      RETURN
      END
