      program   driver
c     ------------------------------------------------------------------
c     Fortran driver for subroutine gqapd (A GRASP for QAP).
c     ------------------------------------------------------------------
c     Parameters, variables and arrays needed to use subroutine gqapd.
c     ------------------------------------------------------------------
      integer   nmax,nsq,in,out
      parameter (nmax=256,nsq=nmax*nmax,in=5,out=6)
      integer   f(nsq),d(nsq),a(nmax),b(nmax)
      integer   srtf(nsq),srtif(nsq),srtd(nsq)
      integer   srtid(nsq),srtc(nsq),srtic(nsq)
      integer   indexd(nsq),indexf(nsq),cost(nsq)
      integer   fdind(nsq),perm(nmax)
      integer   itr,maxitr,bestv,n,look4,seed,n2
      real      alpha,beta
c     ------------------------------------------------------------------
c     Variables used by driver (not needed by gqapd).
c     ------------------------------------------------------------------
      integer   iseed0

c     ------------------------------------------------------------------
c     Read problem data and gqapd parameters.
c     Write out summary of input.
c
c     Read problem dimension (n).
c     ------------------------------------------------------------------
      read (in,*) n
      if (n .gt. nmax) then
	   write(out,10)
10         format(' error: n > nmax  (increase nmax in drive.f)')
	   stop
      else
           n2=n*n
c          -------------------------------------------------------------
c          Read matrices D and F.
c          -------------------------------------------------------------
           call readp(in,n2,d)
           call readp(in,n2,f)
      endif
c     
c     ------------------------------------------------------------------
c     Run GRASP for QAP.
c     ------------------------------------------------------------------
c
c     Variables needed for input:
c       
c         n      : dimension of QAP
c                  integer
c         n2     : n * n
c         maxitr : maximum number of grasp iterations
c                  integer
c         alpha  : grasp construction phase parameter alpha (0.5)
c                  real  
c         beta   : grasp construction phase parameter beta (0.1)
c                  real  
c         look4  : grasp returns permutation if solution with cost
c                  less than or equal to look4 is found (look4 = -1
c                  causes grasp to take maxitr iterations
c                  integer
c         seed   : seed for random number generator (270001)
c                  integer
c
c     Integer arrays needed for input:
c
c         f      : flow matrix (represented as row-by-row array of
c                  dimension NMAX*NMAX)
c         d      : distance matrix (represented as row-by-row array of
c                  dimension NMAX*NMAX)
c         
c     Integer arrays needed for work:
c   
c         a      : dimension NMAX
c         b      : dimension NMAX
c         srtf   : dimension NMAX*NMAX
c         srtif  : dimension NMAX*NMAX
c         srtd   : dimension NMAX*NMAX
c         srtid  : dimension NMAX*NMAX
c         srtc   : dimension NMAX*NMAX
c         srtic  : dimension NMAX*NMAX
c         indexd : dimension NMAX*NMAX
c         indexf : dimension NMAX*NMAX
c         cost   : dimension NMAX*NMAX
c         fdind  : dimension NMAX*NMAX
c         
c     Integer array needed for output:
c   
c         perm   : permutation vector of best solution found 
c                  dimension NMAX
c
c     Integer variables needed for output:
c   
c         bestv  : cost of best assignment found
c         itr    : number of grasp iterations taken  
c  
c     ------------------------------------------------------------------
c     Set grasp parameters.
c     ------------------------------------------------------------------
      maxitr=1000
      alpha=.25
      beta=.5
      look4=-1
      seed=270006
      iseed0=seed

c     ------------------------------------------------------------------
c     Find approximate solution to QAP.
c     ------------------------------------------------------------------
      call gqapd(n,n2,maxitr,alpha,beta,look4,seed,f,d,a,b,
     +           srtf,srtif,srtd,srtid,srtc,srtic,indexd,
     +           indexf,cost,fdind,perm,bestv,itr)

c     ------------------------------------------------------------------
c     Write output summary.
c     ------------------------------------------------------------------
      call outsol(out,n,itr,iseed0,bestv,maxitr,look4,alpha,beta,perm)
c     ------------------------------------------------------------------
      stop
      end

 


      subroutine readp(in,n2,array)
c     ------------------------------------------------------------------
c     readp:  Read in from unit in integer array of size n2
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c          in     - input unit
c          n2     - square of QAP dimension (n * n)
c
c     ------------------------------------------------------------------
      integer in,n2
c     ------------------------------------------------------------------
c     Passed output array:
c
c          array  - integer array of size n2
c
c     ------------------------------------------------------------------
      integer array(n2)
c     ------------------------------------------------------------------
c     Local scalar:
c
c          i      - index
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
      read (in,*) array

      return
      end





      subroutine outsol(out,n,iter,iseed0,opt,niter,look4,
     +                  alpha,beta,opta)
c     ------------------------------------------------------------------
c     Output best solution found and run statistics.
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c          out    - output unit
c          n      - QAP dimension
c          iter   - number of GRASP iterations taken
c          iseed0 - initial random number generator seed
c          opt    - value of best permutation found
c          niter  - maximum number GRASP iterations
c          look4  - value of permutation sought
c          alpha  - GRASP construction phase parameter
c          beta   - GRASP construction phase parameter
c          
c     ------------------------------------------------------------------
      integer out,n,iter,iseed0,opt,niter,look4
      real    alpha,beta
c     ------------------------------------------------------------------
c     Passed input array:
c
c          opta   - best permutation found
c
c     ------------------------------------------------------------------
      integer opta(n)
c     ------------------------------------------------------------------
c     Local scalar:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
      write(out,10)
10    format(
     +' ----------------------------------------------------------')
      write(out,20)
20    format(
     +' G R A S P for Q A P---------------------------------------')
      write(out,30)
30    format('                                        ')

      write(out,40)
40    format(
     +' input-----------------------------------------------------')
      write(out,50) n
50    format('  dimension of qap                   : ',i20)
      write(out,60) alpha
60    format('  construction phase parameter alpha : ',f20.2)
      write(out,70) beta
70    format('  construction phase parameter beta  : ',f20.2)
      write(out,80) niter
80    format('  maximum number of grasp iterations : ',i20)
      write(out,90) look4
90    format('  look4                              : ',i20)
      write(out,100) iseed0
100   format('  initial seed                       : ',i20)
      write(out,30)

      write(out,110)
110   format(
     +' output----------------------------------------------------')
      write(out,120) iter
120   format('  grasp iterations                   : ',i20)
      write(out,140) opt
140   format('  cost of best permutation found     : ',i20)
      write(out,150) (opta(i),i=1,n)
150   format('  best permutation found             : ',
     &       5i4/(37x,': ',5i4))
      write(out,30)
      write(out,10)
      print *,opt
      print *,(opta(i),i=1,n)
      return
      end

