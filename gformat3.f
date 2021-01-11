C     PROGRAM TO CONVERT 32 DETECTOR BT1 FILES INTO SINGLE HISTOGRAM GSAS FILES
C     Original: REM 6/93
C     Modifications: J. Stalick
C     Reorganized: B. Toby
C     Updated 9/07:  JKS
C     
      INTEGER MAX_IN,MAX_OUT
      PARAMETER (MAX_IN=2000)
      PARAMETER (MAX_OUT=99999)
C     input data
      DIMENSION FINT(32,MAX_IN)
      DIMENSION SIG(32,MAX_IN)

C     normalization
      REAL PRESET,PRESET1,PNORM
C     rebinned data
      real WSINT(32,MAX_IN),sumW(32,MAX_IN),TTSTART(32)
      INTEGER NPT(32)
C     overlapped differences
      DIMENSION AVINTDIFF(32),FCORR(32)

      DIMENSION FINTB(32,MAX_IN)
      DIMENSION SIGB(32,MAX_IN)

      DIMENSION COEF(20)
C     pseudo-single detector arrays
      DIMENSION FINTV(MAX_OUT)
      DIMENSION SIGV(MAX_OUT)
      REAL CHEB(30)

      CHARACTER*80 GCARD,LINE
      CHARACTER*72 CARD
      CHARACTER*5 MONO
      CHARACTER*51 TITLE
      CHARACTER*10 TEN
      CHARACTER*4 ID
      CHARACTER*1 IREF(20)
      integer lench
      integer iargc
      logical errorflag,splineflag,interpolate,overwrite,sumflag,
     $     pipeflag,filterflag,absflag,rescale
      integer i,j,nargs,nfiles,ln,lastdet,glast,NTMP
      character*50 arg,filename,outfile,parmfile
      character*20 files(50)
      integer statcheck
C     
C     FIRST SET UP THE INITIAL PARAMETERS
C     
      NTERMR = 5
      NIN = 1
      NOUT = 7
      NINST = 8
      NTMP = 9
      statcheck = 0
      open(unit=9,file='stat.check',status='unknown',access='append')
C     
C     get commandline options--------------------------------------------------
      nargs =  iargc()
      splineflag = .true.
      interpolate = .true.
      errorflag = .false.
      overwrite = .false.
      pipeflag = .false.
      filterflag = .false.
      sumflag = .false.
      absflag = .false.
      rescale = .true.
      lastdet = 32
      nfiles = 0
C     loop over arguments
      do I=1,nargs
         call getarg(I,arg)
         ln = lench(arg) 
C     is this a flag?
         if (arg(1:2) .eq. '-l' .or. arg(1:2) .eq. '-L') then
            read(arg(3:),*,iostat=j) lastdet
            if (j .ne. 0) then
               write(0,*) 'Unable to parse: ',arg
               errorflag = .true.
            elseif (lastdet .le. 0 .or. lastdet .gt. 32) then
               write(0,*) 'Invalid detector number in: ',arg,lastdet
               errorflag = .true.
            endif
         elseif (arg(1:1) .eq. '-') then
            do j=2,ln
               if (arg(j:j) .eq. 'o' .or. arg(j:j) .eq. 'O') then
                  splineflag = .false.
               elseif (arg(j:j) .eq. 'm' .or. arg(j:j) .eq. 'M') then
                  interpolate = .false.
               elseif (arg(j:j) .eq. 'c' .or. arg(j:j) .eq. 'C') then
                  overwrite = .true.
               elseif (arg(j:j) .eq. 's' .or. arg(j:j) .eq. 'S') then
                  sumflag = .true.
               elseif (arg(j:j) .eq. 'f' .or. arg(j:j) .eq. 'F') then
                  filterflag = .true.
		  pipeflag = .true.
               elseif (arg(j:j) .eq. 'a' .or. arg(j:j) .eq. 'A') then
                  absflag = .true.
               elseif (arg(j:j) .eq. 'r' .or. arg(j:j) .eq. 'R') then
                  rescale = .false.
               else
                  write(0,'(/,5x,2a,/)')
     $                 'Invalid command line flag: -',arg(j:j)
                  nargs = 0
                  goto 100
               endif
            enddo
         else
C     remove the .raw, if present
            if (arg(ln-3:ln) .eq. '.raw') then 
               ln = ln - 4
            endif
C     is this a valid file?
            open (unit=1,file=arg(1:ln)//'.raw',status='old',
     $           iostat=j)
            if (j .ne. 0) then
               write(0,'(/,5x,2a)') 'Error opening file ',
     $              arg(1:ln)//'.raw'
               errorflag = .true.
            else
               nfiles = nfiles + 1
               if (nfiles .gt. 50) then
                  write(0,*)
     $                 'Sorry, this program can only process 50 files'
                  stop
               endif
               files(nfiles) = arg(:ln)
            endif
            close(unit=1)
         endif
C     end of loop over arguments
      enddo

      if (nfiles .eq. 1 .and. sumflag) then
         write (0,*) 'Only 1 file to process, ignoring -s flag'
         sumflag = .false.
      endif
      if (nfiles .ne. 0 .and. filterflag) then
         write (0,*) 'You cannot specify -f and a input file'
         filterflag = .false.
         nfiles = 0
      elseif (filterflag) then
         nfiles = 1
      endif
 100  if (nargs .le. 0 .or. errorflag .or. nfiles .le. 0) then
         write(0,'(/,1x,a)')
     $        'Transform BT1 data to pseudo-single'//
     $        ' detector data for GSAS'
         write(0,'(/,1x,a)') 'Usage:',
     $        '  gformat [-Lxx] [-coms] file1[.raw] [file2[.raw]] ...',
     $        '  gformat -f [-Lxx] [-om] < file1.raw > file2.gsas'
         write(0,'(9x,3a)') '-c: ',
     $        ' overwrite (clobber) existing output files'
         write(0,'(9x,3a)') '-m:  move points to closest bin',
     $        '[default: interpolate]'
         write(0,'(9x,2a)') '-s: ',
     $        'sum multiple files '//
     $        '[default: process each file independently'
         write(0,'(9x,3a)') '-o: ',
     $        ' omit background spline',
     $        ' [default: add offsets to spline backgd]'
         write(0,'(9x,2a)')
     $        '-f: run as filter: read from stdin and write to stdout'
         write(0,'(9x,2a)') '-Lxx or -lxx: ',
     $        'drop detectors xx-32 [default: use all]'
         write(0,'(9x,2a,/)') '-a: ',
     $        'apply absorption correction from file abs.corr'
         write(0,'(9x,2a,/)') '-r: ',
     $        'prevent rescaling to ~counts'
         stop
      endif
C     done getting commandline options-----------------------------------------
      WRITE (0,10)
 10   FORMAT ('Transformation of 32 det BT1 data to GSAS format',/)

C     status messages
      if (splineflag) then
         write (0,*)
     $        'Offsets will be added to data to spline backgrounds'
      else
         write (0,*)
     $        'No spline correction will be applied to backgrounds'
      endif
      if (interpolate) then
         write (0,*) 'Data will be interpolated'
      else
         write (0,*) 'Data will not be interpolated'
      endif
      if (overwrite) then
         write (0,*) 'Will overwrite existing output files!'
      endif

C-----------------------------------------------------------------------
C     loop over files       ii is the file number
C-----------------------------------------------------------------------
      do ii=1,nfiles
C     open the input file
	 if (.not. filterflag) then
            filename = files(ii)(:lench(files(ii))) // '.raw'
            OPEN (UNIT=NIN,  FILE= filename, STATUS= 'OLD',
     $           iostat=j)
            if (j .ne. 0) then
               write(0,'(/,5x,2a)') 'Error opening file ',
     $              files(ii)(:lench(files(ii))) // '.raw'
               stop
            endif
            WRITE (0,'(a)') ' Reading input from '//filename
C     put the filename into the statcheck file, always
            write (9,*) 'file ',files(ii)(:lench(files(ii))) // '.raw'
         else
            filename = 'stdin'
            WRITE (0,'(a)') ' Reading input from stdin'
            nin = 5
         endif
         if (ii .eq. 1 .or. .not. sumflag) then
C     zero the arrays
            DO N=1,MAX_IN
               DO J=1,32
                  WSINT(J,N) = 0
                  sumW(J,N) = 0
               ENDDO
            ENDDO
            DO J=1,32
               npt(j) = 0
            ENDDO
         endif
C     open output file
         if ((ii .eq. 1 .or. .not. sumflag) .and. (.not. pipeflag)) then
            errorflag = .false.
            outfile = files(ii)(:lench(files(ii))) // '.gsas'
            parmfile = files(ii)(:lench(files(ii))) // '.inst'
            OPEN (UNIT=NOUT, FILE= outfile, STATUS= 'NEW', iostat=j)
            if (overwrite .and. j .ne. 0) THEN
               OPEN (UNIT=NOUT, FILE= outfile, STATUS= 'OLD', iostat=j)
               if (j .ne. 0) then
                  WRITE(0,'(3a)') '  ***ERROR***  Unable to overwrite ',
     $                 outfile(:lench(outfile)),
     $                 ' Protection problem?'
                  errorflag = .true.
               ENDIF
            elseif (j .ne. 0) then
               WRITE(0,'(3a)') '  ***ERROR***  Unable to create ',
     $              outfile(:lench(outfile)),
     $              ' Does it already exist?'
               errorflag = .true.
            endif
            if (.not. errorflag) write(0,*)
     $           'Creating ',outfile(:lench(outfile)) 
            if (errorflag) stop
         elseif (ii .eq. 1 .or. .not. sumflag) then
            WRITE (0,'(a)') ' Sending output to stdout'
            nout = 6
         endif
C     
C     
C     READ IN FIRST TWO LINES OF INPUT FILE
C     
C     
         READ (nin,'(a50,f8.0,2i3,1x,a5,g14.1)')
     $        TITLE,WL1,icol1,icol2,mono,preset
         if (ii .eq. 1 .or. .not. sumflag) then
            PRESET1 = PRESET
         ENDIF
         IF (sumflag .and. (PRESET .EQ. 0 .or. PRESET1 .eq. 0)) THEN
            PNORM = 1.0
            WRITE(0,*) 'No prefactor in ',
     $           files(ii)(:lench(files(ii))) // '.raw',
     $           '. Rerun proprep to sum.'
         ELSEIF (sumflag) THEN
            PNORM = PRESET1 /  PRESET
         ELSE
            PNORM = 1.0
         ENDIF
         write (0,*) 'PNORM=',pnorm

         READ (NIN,'(A72)') CARD
         ID = CARD(1:4)
!     write(0,*) 'colls = ',icol1,icol2
!     write(0,*) 'mono  = ',mono
C     
C-----OUTPUT HEADER #1
C     
C     create the title line
         if (ii .eq. 1 .or. .not. sumflag) then
            GCARD = ' '
            GCARD(1:51) = TITLE
            glast = 80
         endif
C     add the name of each file into the header
         ln = lench(files(ii))
         if (glast .gt. 20) gcard(glast-ln:glast) = files(ii)(:ln)
!     write (*,*) glast
!     write (*,*) gcard,ln,files(ii)(:ln)

         glast = glast - ln - 1
C     
C     
C     READ IN THE INTENSITIES, AND CALCULATE THE CORRECT 2THETA 
C     VALUES ASSOCIATED WITH THEM
C     
C     
         NUMDET = 0
C-----------------------------------------------------------------------
C     loop over detectors      J is the detector number
C-----------------------------------------------------------------------
         do j=1,lastdet
C     
C     READ BANK HEADER
            READ (NIN,'(F8.0,F8.0,F8.0,2F8.4,2I5)',err=412)
     $           TMIN, FSTEP, TMAX, FMULT, ZERO,IDET,IPT
            TMIN = TMIN/100.
            TMAX = TMAX/100.
            FSTEP = FSTEP/100.
            ZERO = ZERO/100.
            NUMDET = NUMDET + 1
            if (ipt .gt. MAX_IN) then
               write(0,*) 'Number of input points greater than MAX_IN',
     $              ipt,max_in
               stop
            endif
            
            READ (NIN,'(12F6.0)') (FINT(J,I),I=1,ipt)

C     
C-----CORRECT INTENSITIES and esd's USING DETECTOR SCALES (FMULT)
C     
            DO I = 1,IPT
               if (fint(j,i) .gt. 0)
     $              sig(J,i) = fmult * sqrt(fint(J,i)) * PNORM
               FINT(J,I) = FINT(J,I) * FMULT * PNORM
               IF (fint(j,i) .eq. 0) sig(J,i) = 1.
            enddo

C     1st file, 1st detector
            IF (ii .eq. 1 .and. j .eq. 1 .or. .not. sumflag) then
               IF (INTERPOLATE) THEN
C     GSTEP is output step size (output=input for interpolation, doubled otherwise)
C     this could be a user controlled option someday
                  GSTEP = FSTEP
               ELSE
!     GSTEP = 2.*FSTEP
                  GSTEP = FSTEP
               ENDIF
            ENDIF

C     1st file, 1st point (all detectors) -- save starting angle
            IF (ii .eq. 1 .or. .not. sumflag) then
C     correct 2theta
               FTHETA = TMIN - ZERO
               IF (INTERPOLATE) THEN
                  REMAIN = FTHETA - INT(FTHETA/GSTEP)*GSTEP
                  N = 1
C     ignore a difference of less than 0.01 centidegree (plus 0.00001 for roundoff)
                  IF (ABS(REMAIN) .LT. 0.00011) N = 0
C     compute the 2theta that is an even multiple of the output stepsize and 
C     larger than the starting angle
                  TTSTART(J)  = (INT(FTHETA/GSTEP)+N)*GSTEP
               ELSE
C     compute a 2theta that is an even multiple of the output stepsize)
                  TTSTART(J) = (NINT(FTHETA/GSTEP))*GSTEP
               ENDIF
!     write (*,*) j,' 2theta, bin start', FTHETA, TTSTART(J)
            ENDIF
C     
C-------CORRECT INTENSITIES FOR BIN--Assuming straight lines 
C     successive points
C     
C     
            IF (INTERPOLATE) IPT = IPT - 1
            DO I = 1,IPT
C     correct 2theta
               FTHETA = TMIN - ZERO + (I-1)*FSTEP
               IF (INTERPOLATE) THEN
C     
C------CALCULATE THE RIGHT THETA VALUE FOR GSAS BIN
C     THIS ROUNDS UP TO THE NEAREST MULTIPLE OF GSTEP
                  REMAIN = FTHETA - INT(FTHETA/GSTEP)*GSTEP
                  N = 1
C     ignore a difference of less than 0.01 centidegree (plus 0.00001 for roundoff)
                  IF (ABS(REMAIN) .LT. 0.00011) N = 0
C     compute a 2theta that is an even multiple of the output stepsize)
                  FTHETAC = (INT(FTHETA/GSTEP)+N)*GSTEP
                  RATIO = (FTHETAC - FTHETA)/FSTEP
                  FINTERP = FINT(J,I) +
     $                 (FINT(J,I+1)-FINT(J,I)) * RATIO
                  sigintepb = sig(J,i) + 
     $                 (sig(J,i+1)-sig(J,i)) * RATIO
				  c1 = (1.-1./RATIO) * (1.-1./RATIO)
				  c2 = 1./RATIO/RATIO
				  sigj2 = sig(J,i) * sig(J,i)
				  sigj2b = sig(J,i+1) * sig(J,i+1)
                  sigintep = sqrt( c1 * sigj2 + c2 * sigj2b)
                  
                  WRITE (0,*) 'sig is times',sigintep
                  WRITE (0,*) 'sigb is times',sigintepb

C     deal with negative intensities (ignore them)
                  IF (FINT(J,I) .lt. 0 .or. FINT(J,I+1) .lt. 0)
     $                 sigintep = 0.0
               ELSE
C     no interpolation
C------CALCULATE THE CLOSEST GSAS BIN
                  FTHETAC = (NINT(FTHETA/GSTEP))*GSTEP
                  FINTERP = FINT(J,I)
                  sigintep = sig(J,i)
C     deal with negative intensities (ignore them)
                  IF (FINT(J,I) .lt. 0) sigintep = 0.0
               ENDIF
C     compute the point number for this 2theta
               N = NINT((FTHETAC - TTSTART(J))/GSTEP) + 1
               NPT(J) = MAX(NPT(J),N)
               if (sigintep .gt. 0) then
                  s2 = 1./sigintep**2
C     Stat check#1 -- one detector to previous measurements of the same detector
                  IF (sumW(J,N) .gt. 0) then
                     diff =  FINTERP-WSINT(J,N)/sumW(J,N)
                     sigdif = sqrt (1./sumW(J,N) + 1./s2)
                     IF (abs(diff) .gt. 4.*sigdif) then
                        statcheck = statcheck+1
                        write (9,*)
     $     '  Warning: file not consistent with previous'
     $     ,' by ',nint(10*abs(diff)/sigdif)/10.,' sigma'
                        write (9,*) 'detector ',J,' angle',FTHETA
                        write (9,*) 'Intensities',
     $                       WSINT(J,N)/sumW(J,N),FINTERP
                        write (9,*) 'ESD''s      ',
     $                       sqrt(1./sumW(J,N)),sigintep 
                     endif
                  ENDIF
C     weighted sum of intensities
                  WSINT(J,N) = WSINT(J,N) + FINTERP*s2
C     sum of weights
                  sumW(J,N) = sumW(J,N) + s2
               ENDIF
            ENDDO
C     End loop over all detectors
         ENDDO
 412     CLOSE(NIN)

C     process and write data unless there are more files to be summed
         IF (II .eq. nfiles .or. .not. SUMFLAG) then
C     
C-----Correct for different backgrounds in different
C-----detectors. Use detector #20 as a reference baseline.
C     
C-----Find the points in the previous detector that overlap the current
C-----Start with the second detector
            DO J = 1,lastdet-1
               OVERSUM1 = 0.0
               OVERWT1 = 0.0
               OVERSUM2 = 0.0
               OVERWT2 = 0.0
C     compute the point number for the first overlapped point
               N = NINT((TTSTART(J+1) - TTSTART(J))/GSTEP) + 1
               I1 = 0
C     write(*,*) NPT(J)-N+1,' Points',(NPT(J)-N+1)*GSTEP,' deg.'
               DO I=N,NPT(J)
                  I1 = I1 + 1
                  OVERSUM1 = OVERSUM1 + WSINT(J,I)
                  OVERSUM2 = OVERSUM2 + WSINT(J+1,I1)
                  OVERWT1 = OVERWT1 + sumW(J,I)
                  OVERWT2 = OVERWT2 + sumW(J+1,I1)
               ENDDO
               OVERSUM1 =OVERSUM1/OVERWT1
               OVERWT1 = 1./SQRT(OVERWT1)
               OVERSUM2 =OVERSUM2/OVERWT2
               OVERWT2 = 1./SQRT(OVERWT2)
!     write(*,'(2(a,i3,f8.2,a,f7.2))') 'overlap ',J,
!     $              OVERSUM1,' +-',OVERWT1,' compared to ',j+1,
!     $              OVERSUM2,' +-',OVERWT2
C     Is this difference statistically significant? 
C     It must be greater than lets say 2 sigma
               IF (ABS(OVERSUM1-OVERSUM2) .gt.
     $              2*SQRT(OVERWT1**2 + OVERWT2**2)) THEN
                  AVINTDIFF(J) = OVERSUM1 - OVERSUM2 
!     write(*,*) 'correction=',
!     $   OVERSUM1-OVERSUM2,' +-',SQRT(OVERWT1**2 + OVERWT2**2),
!     $                 ' is meaningful'
               ELSE
                  AVINTDIFF(J) = 0.0
!     write(*,*) 'correction=',
!     $   OVERSUM1-OVERSUM2,' +-',SQRT(OVERWT1**2 + OVERWT2**2),
!     $                 ' not meaningful, not used'
               ENDIF
            ENDDO
C     
C---- Calculate the actual corrections that need 
C---- to be applied to each detector.
C---- Detector 20 = baseline
C     
C     note that Fcorr will be added to the appropriate detector
            do j=1,32
               fcorr(j) = 0.
            enddo

            IF (SPLINEFLAG) THEN
c     
c---- Then first 19 detectors
c     
               do  i = 19,1,-1
                  fcorr(i) = fcorr(i+1) - AVintdiff(i)
               enddo
C     
C---- then remaining detectors
C     
               do i = 21,lastdet
                  fcorr(i) = fcorr(i-1) + AVintdiff(i-1)
               enddo
C     tell the user what we will do
               write(0,'(2a)') 'The following background spline ',
     $              'corrections will be added to each detector'
               write(0,'(6(i4,a,f8.2))') (j,':',fcorr(j),j=1,lastdet)
            ENDIF

C     
C     
C---- Normalize the intensities and apply these corrections to the detectors
C     
C     
            do  i = 1,lastdet
               do  j = 1,NPT(I)
                  if (sumW(I,J) .gt. 0) THEN
                     fintb(I,j) = wsint(I,j)/sumW(I,J) + fcorr(i)
                     sigb(I,J) = 1./SQRT(sumW(I,J))
                  ELSE
                     fintb(I,j) = 0.0
                     sigb(I,J) = 0.0
                  ENDIF
C     write (13,445) TTSTART(I)+(J-1.)*GSTEP,fintb(i,j)
C     445       format(f8.3,f8.1)
               enddo
            enddo
            

C     compute range of data and zero array
            TSTART = TTSTART(1)
            TEND = TTSTART(LASTDET) + (NPT(LASTDET)-1.)*GSTEP
C     write (*,*) TTSTART(1),TTSTART(LASTDET),NPT(LASTDET)
            NPNTS = 1 + NINT( (TEND - TSTART)/GSTEP )
            IF (NPNTS .gt. MAX_OUT) THEN
               write(0,*) 'Warning: # points > MAX_OUT',
     $              npnts,MAX_OUT
               write(0,*) 'Not all data will be written'
               NPNTS = MAX_OUT
            ENDIF
            DO I=1,NPNTS
               FINTV(I) = 0.
               SIGV(I) = 0.
            ENDDO
C     Now put into pseudo-single detector array
            DO J=1,LASTDET
               DO I=1,NPT(J)
C     compute point number in master array
                  N = I + NINT((TTSTART(J) - TTSTART(1))/GSTEP)
                  IF (sigb(j,i) .gt. 0 .and. N .le. MAX_OUT) THEN
                     s2 = 1./sigb(j,i)**2 
C     Stat check#2 -- current detector to previous one
                     IF (SIGV(N) .gt. 0) then
                        diff =   FINTB(j,i) - FINTV(N)/SIGV(N)
                        sigdif = sqrt (1./SIGV(N) + 1./s2)
                        IF (abs(diff) .gt. 4.*sigdif) then
                           statcheck = statcheck+1
                           write (9,*)
     $    '  Warning: detector not consistent with previous'
     $    ,' by ',nint(10*abs(diff)/sigdif)/10.,' sigma'
                           write (9,*) 'detector ',J,' angle',
     $                          GSTEP*(N-1) + TTSTART(1)
                           write (9,*) 'Intensities',
     $                          FINTV(N)/SIGV(N),FINTB(j,i)
                           write (9,*) 'ESD''s      ',
     $                          sqrt(1./SIGV(N)),sigb(j,i)
                        endif
                     ENDIF
                     FINTV(N) = FINTV(N) + FINTB(j,i) * s2
                     SIGV(N) = SIGV(N) + s2
                  ENDIF
               ENDDO
            ENDDO
C     normalize
            ICOUNT = 0
            DO N=1,NPNTS
               IF ( SIGV(N) .GT. 0) THEN
                  FINTV(N) = FINTV(N) / SIGV(N)
                  SIGV(N) = 1./SQRT(SIGV(N))
               ELSE
                  ICOUNT = ICOUNT+1
               ENDIF
            ENDDO
            IF (ICOUNT .gt. 0) WRITE(0,*) 'Warning ',ICOUNT,
     $           ' Points were not defined'
            IF (ABSflag) THEN
C     absorption correction
               WRITE (0,*) 'Applying Absorption Correction'
               OPEN(unit=99,file='abs.corr',status='old')
C     file abs.corr contains a one line header
C     and then the number of Chebyshev terms, followed by the terms
               read(99,'(A)') title
               read(99,*) NTERM,(CHEB(I),I=1,NTERM)
               close(unit=99)
               WRITE (0,*) '  from File abs.corr', title
               WRITE (0,*) NTERM,'coef: ',(CHEB(I),I=1,NTERM)
               absmin = 1e9
               absmax = 0
               DO N=1,NPNTS
C     compute 2theta in degrees
                  TTHETA = TSTART + (N-1.)*GSTEP
C     change range to +-1
                  T = -1. + TTHETA/90.
                  DY1 = 1.
                  DY2 = T
                  abscor = CHEB(1) + T*CHEB(2)
                  DO I=3,NTERM
                     DYI = 2. * T * DY2-DY1
                     DY1 = DY2
                     DY2 = DYI
                     abscor = abscor + CHEB(I)*DYI
                  ENDDO
                  absmin = min(absmin, abscor) 
                  absmax = max(absmax, abscor) 
                  FINTV(N) = FINTV(N) * ABSCOR
                  SIGV(N) = SIGV(N) * ABSCOR
               ENDDO
               WRITE (0,*) 'Correction ranged from ',absmin,
     $              ' to ', absmax 
            ENDIF

            IF (RESCALE) THEN
               ynorm = 0
               icount = 0
               DO N=1,NPNTS
                  IF ( SIGV(N) .GT. 0) THEN
                     ynorm = ynorm + FINTV(N) / (SIGV(N)*SIGV(N))
                     ICOUNT = ICOUNT+1
                  ENDIF
               ENDDO
               IF (ICOUNT .GT. 0) THEN
                  YNORM = YNORM / ICOUNT
               ELSE
                  YNORM = 1
               ENDIF
               WRITE (0,*) 'Rescaling factor is times',ynorm
               DO N=1,NPNTS
                  IF ( SIGV(N) .GT. 0) THEN
                     FINTV(N) = ynorm * FINTV(N)
                     SIGV(N)  = ynorm * SIGV(N)
                  ENDIF
               ENDDO
            ENDIF
C     
C---- OUTPUT TO GSAS FORMAT FILE (NOUT)
C---- NPNTS = NO OF GSAS POINTS
C---- CALCULATE NO OF GSAS RECORDS
C     
            IREC = INT((NPNTS-1)/5.)+1

            sstart = INT(100*TSTART)
            sstep = INT(100*GSTEP)
            WRITE(NOUT,'(3A)') GCARD
            if (.not. pipeflag) then
               WRITE(LINE,'(A26,A)') 'Instrument parameter file:',
     $              parmfile(:lench(parmfile))
               WRITE(NOUT,'(3A)') LINE
            endif
            WRITE(LINE,500) NPNTS,IREC,sstart,sstep
 500        FORMAT('BANK  1 ',2i8,' CONST ',2f10.2,' 0 0 ESD',21x)
            WRITE(NOUT,'(3A)') LINE
            N = 1
            DO 530 I = 1,IREC
               WRITE (LINE,520) (FINTV(J),sigv(j),J=N,N+4)
 520           FORMAT (5(F8.0,f8.0))
               WRITE(NOUT,'(3A)') LINE
               N = N + 5
 530        CONTINUE 
C     
C     
            if (.not. pipeflag) then
               OPEN (UNIT=NINST, FILE= parmfile, STATUS= 'unknown',
     $              err=919)
               write(0,*) 'Creating ',parmfile(:lench(parmfile)) 
C-----OUTPUT INSTRUMENT PARAMETER FILE (NINST)
C     
C     OUTPUT FIRST 3 LINES TO INSTRUMENT PARAMETER FILE
C     
               TEN = '1234567890'
               WRITE (LINE,535) TEN,TEN,TEN,TEN,TEN,TEN
 535           FORMAT (12X,6A10)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,540) 1
 540           FORMAT ('INS   BANK  ',I5)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,550)
 550           FORMAT ('INS   HTYPE   PNCR')
               WRITE(NINST,'(3A)') LINE
C     
C     
C     
C-----MISCELANEOUS SET UP PARAMETERS
               DO  I = 1,20
                  IREF(I) = 'N'
               ENDDO
               DO I = 5,20
                  COEF(I) = 0.0
               ENDDO
               IF (mono .eq. 'GE311') THEN
C     Ge311 parameters
                  if (icol1 .eq. 7) then
                     COEF(1) = 246.7 
                     COEF(2) = -298.6
                     COEF(3) = 157.6
                     COEF(4) = 4.0305
                  else
                     COEF(1) = 398.5
                     COEF(2) = -343.2
                     COEF(3) = 163.0
                     COEF(4) = 4.0305
                  endif
               ELSEIF  (mono .eq. 'GE733') THEN
C     Ge733 parameters
                  if (icol1 .eq. 7) then
                     COEF(1) = 36.0 
                     COEF(2) = -100.7 
                     COEF(3) = 106.0
                     COEF(4) = 8.0
                  else
                     COEF(1) = 51.3
                     COEF(2) = -109.0
                     COEF(3) =  105.6
                     COEF(4) = 8.0
                  endif
               ELSE
C     Cu311 parameters
                  if (icol1 .eq. 7) then
                     COEF(1) = 181.2
                     COEF(2) = -299.3
                     COEF(3) = 186.3
                     COEF(4) = 5.0
                  else
                     COEF(1) = 239.7
                     COEF(2) = -298.2
                     COEF(3) = 180.8
                     COEF(4) = 5.0
                  endif
               ENDIF
C     
               WRITE (LINE,555) WL1,0.0,ZERO,0
 555           FORMAT ('INS  1 ICONS',3F10.5,I10,10X)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,560) ID
 560           FORMAT ('INS  1I HEAD',2X,A4,2X,60X)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,565) 0,TSTART,TMAX,1
 565           FORMAT ('INS  1I ITYP',I5,2F10.4,I10)
               WRITE(NINST,'(3A)') LINE

               IPTYP=1
               IDAMP = 0
               WRITE (LINE,575) IPTYP,IPTYP,6,0.005,IDAMP,
     $              (IREF(I),I=1,6)
 575           FORMAT ('INS  1PRCF',I1,1X,2I5,F10.5,4X,I1,20A1)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,580) IPTYP,1,(COEF(I),I=1,4)
 580           FORMAT ('INS  1PRCF',2I1,4E15.6)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,580) IPTYP,2,(COEF(I),I=5,6)
               WRITE(NINST,'(3A)') LINE
               IPTYP=2
               COEF(7) = COEF(4)
               COEF(4) = 0.0
               WRITE (LINE,575) IPTYP,IPTYP,12,0.005,IDAMP,
     $              (IREF(I),I=1,12)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,580) IPTYP,1,(COEF(I),I=1,4)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,580) IPTYP,2,(COEF(I),I=5,8)
               WRITE(NINST,'(3A)') LINE

               WRITE (LINE,580) IPTYP,3,(COEF(I),I=9,12)
               WRITE(NINST,'(3A)') LINE
               IPTYP=3
               COEF(7) = 0.04
               COEF(8) = 0.03
               WRITE (LINE,575) IPTYP,IPTYP,19,0.005,IDAMP,
     $              (IREF(I),I=1,19)
               WRITE(NINST,'(3A)') LINE
               
               WRITE (LINE,580) IPTYP,1,(COEF(I),I=1,4)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,580) IPTYP,2,(COEF(I),I=5,8)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,580) IPTYP,3,(COEF(I),I=9,12)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,580) IPTYP,4,(COEF(I),I=13,16)
               WRITE(NINST,'(3A)') LINE
               WRITE (LINE,580) IPTYP,5,(COEF(I),I=17,19)
               WRITE(NINST,'(3A)') LINE
               CLOSE(ninst)
            endif
            if (.not. pipeflag) then
               CLOSE(nout)
            ENDIF

C     end of if block to write data
 919        CONTINUE
         ENDIF
C     End loop over all files
      ENDDO
      if (statcheck .gt. 0) then
         write(0,*) 'Note: ',statcheck,
     $        ' inconsistency warnings were written to file stat.check'
      endif
      close(9)
      STOP 'OK'
      END
C---------------------------------------------------------------------------
C---------------------------------------------------------------------------
c     Function LENCH
c     
c     This function takes a character string and finds out how long the
c     "actual" string is (i.e. not including padded blanks on the right).
c     
C---------------------------------------------------------------------------
      integer function lench(string)
      character*(*) string
      character*1 nul
      data nul/"\0"/
      lench=0
      if (string.eq.' '.or.string(1:1).eq.nul) return
      do lench=len(string),1,-1
         if (string(lench:lench).ne.' '.and.string(lench:lench).ne.nul)
     $        return
      enddo
      lench=0
      return
      end
C---------------------------------------------------------------------------
C     Convert strings to lower case characters
      SUBROUTINE DOWNCASE(STRING)
      CHARACTER*(*)   STRING    !String to convert
      NCH = LEN(STRING)
      DO ICH=1,NCH
         JCH = ICHAR(STRING(ICH:ICH))
         IF ( JCH.GE.ICHAR('A') .AND. JCH.LE.ICHAR('Z') )
     $        STRING(ICH:ICH) = CHAR(JCH-ICHAR('A')+ICHAR('a'))
      END DO
      RETURN
      END
C---------------------------------------------------------------------------
C     Convert strings to upper case characters
      SUBROUTINE UPCASE(STRING)
      CHARACTER*(*)   STRING    !String to convert
      NCH = LEN(STRING)
      DO ICH=1,NCH
         JCH = ICHAR(STRING(ICH:ICH))
         IF ( JCH.GE.ICHAR('a') .AND. JCH.LE.ICHAR('z') )
     $        STRING(ICH:ICH) = CHAR(JCH-ICHAR('a')+ICHAR('A'))
      END DO
      RETURN
      END
