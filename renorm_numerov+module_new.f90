      MODULE TABS
      IMPLICIT NONE
      COMPLEX*16,ALLOCATABLE,DIMENSION(:,:,:) :: VV,CVV
      COMPLEX*16,ALLOCATABLE,DIMENSION(:) :: integrale_3J 
      COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: POTLM,C,C_inv,C_croix,VV_bis,IDEN,VV_real_bis 
      COMPLEX*16 :: AI,DETC
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: VV_real,CVV_real
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RR,FAC,AC,POTlm_realbis
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: AAA,BBB,CCC
      DOUBLE PRECISION :: RK,RSTART,REND,H,PI,DET,R0,enj,ZEFF
      INTEGER,ALLOCATABLE,DIMENSION(:) :: LPOT,MPOT,L,M
      INTEGER NS,NBR,LMAX,MMAX,NR,NRINI,LPOTMAX,MPOTMAX,IDFA
      PARAMETER(IDFA=300)
      PARAMETER(PI=3.141592654D0)
      END MODULE TABS
      

      PROGRAM DIABATS                                             
      USE TABS
      IMPLICIT NONE
      COMPLEX*16 potentiel,u,V_sum
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: XN,YN,ZN,ZTN
      DOUBLE PRECISION :: PASR,POTLMR,POTLMI,ExR,ExI,RMIN,RMAX,ENJ_W,IP,IP_EV,ENJ_EV,RMAXCHOIX
      INTEGER :: NS_pot,I,J,K,L1,M1,IR
      INTEGER :: M3,N,KK,N_ENJ
      CHARACTER(50) lstring,mstring
      CHARACTER(70) instruction

      AI=dcmplx(0.d0,1.d0)
      
      ALLOCATE(FAC(IDFA),AC(IDFA))
      CALL LOGFAC

	write(6,*) 'coucou test2'
      WRITE(6,7000)                                                
 7000 FORMAT(' ***** PROGRAM RENORM NUMEROV METHOD ***** ')     
      OPEN(UNIT=55,FILE='renorm_numerov.dat',STATUS='old')
      OPEN(UNIT=65,FILE='../RUN/tableau_ener.dat',STATUS='old')

! define wave number for all channels 
      READ(65,*) ENJ_EV       !ENJ_EV est l'energie d'ionisation en ev   
      READ(55,*) IP		! IP en ua
      ENJ_W=ENJ_EV/27.2d0    !ENJ_W est l'energie d'ionisation en ua
      ENJ=ENJ_W-IP           !ENJ est l'energie du photoelectron en ua
      RK=DSQRT(2.D0*ENJ)
      IF(ENJ.le.0.d0)then
      write(6,*) 'IP=',IP,'ENJ_W=',ENJ_W,'pb IP'
      stop 222
      endif

! grid of internuclear distances at which calculations are done         
      READ(55,*) H,RSTART,REND                                      
      WRITE(6,460) H,RSTART,REND                                  
 460  FORMAT(/' H, RSTART, REND ', 3(G12.5,2X))                   

! Lecture des parametres  
      READ(55,*) LMAX,MMAX            !etats pour les fct du continu 
      READ(55,*) LPOTMAX,MPOTMAX      !etats pour le potentiel
      READ(55,*) NRINI                !nb de r que l'on veut pour le potentiel
      READ(55,*) ZEFF                 !charge asymptotique de l'ion

! definition de la grille en R qui sera effectivement utilisee dans l'integration
      N=(REND-RSTART)/H                                             
      H=(REND-RSTART)/DFLOAT(N)
      NR=N-1
      WRITE(6,1000) RSTART,REND,H,NR
 1000 FORMAT(' RSTART=',F10.5,' REND=',F10.5,' H=',F10.5,' NR=',I9)

! Tableau des l et m pour chaque etat (du continu)
      K=0
      DO L1=1,LMAX+1
      DO M1=1,MIN(2*L1-1,2*MMAX+1)
      K=K+1
      ENDDO
      ENDDO
      NS=K
      ALLOCATE(L(NS),M(NS))
      K=0
      DO L1=1,LMAX+1
      DO M1=1,MIN(2*L1-1,2*MMAX+1)
      K=K+1
      L(K)=L1-1
      M(K)=M1-1-MIN(L1-1,MMAX)
      ENDDO
      ENDDO
      WRITE(6,*) 'nombre d"états pour le continu =',NS

! Definition de la grille en r pour le potentiel
      call system('wc -l V_lm_DFT_B3LYP/potentiel_l_0_m_0.dat > nombre_de_r')
      open(8,file='nombre_de_r',status='old')
      read(8,*) nbr
      close(8)

! Ecriture des parametres 
      WRITE(6,*) 'nbr_pot=',NBR
      WRITE(6,*) 'Lmax=',LMAX
      WRITE(6,*) 'Mmax=',MMAX
      WRITE(6,*) 'L_pot_max=',LPOTMAX
      WRITE(6,*) 'M_pot_max=',MPOTMAX


! imax est le nombre d'etats pour le potentiel
      I=0
      DO J=1,LPOTMAX+1
      DO M3=1,MIN(2*J-1,2*MPOTMAX+1)
      I=I+1
      ENDDO
      ENDDO
      NS_pot=I
      WRITE(6,*) 'nombre d"états pour le potentiel =',NS_pot

      ALLOCATE(LPOT(NS_pot),MPOT(NS_pot),POTlm(NRINI,NS_pot))
      ALLOCATE(RR(NRINI),AAA(NRINI),BBB(NRINI),CCC(NRINI))
      ALLOCATE(VV(NRINI,NS,NS),CVV(NRINI,NS,NS))

      I=0
      DO J=1,LPOTMAX+1
      DO M3=1,MIN(2*J-1,2*MPOTMAX+1)
      I=I+1
      LPOT(i)=J-1
      MPOT(i)=M3-1-MIN(j-1,MPOTMAX)
      ENDDO
      ENDDO

! Definition de la matrice de passage aux harmoniques reelles
! les harmoniques reelles sont rangees en colonnes
        allocate(C(NS,NS),C_inv(NS,NS))
        do i=1,NS
        do j=1,NS
! Definition de la matrice par blocs 
        if(l(i).eq.l(j))then
        if((m(i).eq.m(j)).and.m(i).eq.0) then
        C(j,i)=1.d0
        elseif((m(i).eq.m(j)).and.(m(i).lt.0))then
        C(j,i)=(-1.d0)**m(i)*AI/dsqrt(2.d0)
        elseif((m(i).eq.m(j)).and.(m(i).gt.0))then
        C(j,i)=1.D0/dsqrt(2.d0)
        elseif((m(i)+m(j)).eq.0.and.(m(i).lt.0))then
        C(j,i)=-AI/dsqrt(2.d0)
        elseif((m(i)+m(j)).eq.0.and.(m(i).gt.0))then
        C(j,i)=(-1.d0)**m(i)/dsqrt(2.d0)
        else
        C(j,i)=0.d0
        endif
        else
        C(j,i)=0.d0
        endif
        enddo
        enddo
        OPEN(UNIT=999,FILE='MATRICE_C',STATUS='new')
        WRITE(999,*) LMAX,MMAX,NS
        DO I=1,NS
        DO J=1,NS
        WRITE(999,*) DREAL(C(J,I)),DIMAG(C(J,I))
        ENDDO
        ENDDO
        

! definition de la matrice inverse (transposee complexe)
	allocate(IDEN(NS,NS),C_croix(NS,NS))
	do i=1,NS
        do j=1,NS
	C_inv(i,j)=dconjg(C(j,i))
	enddo
	enddo
	IDEN=MATMUL(C_inv,C)
	do i=1,NS
        do j=1,NS
!	write(99,*) i,j,IDEN(i,j)    ! mettre si iden=1 ok sinon stop
	enddo
	enddo

! On lit le potentiel Vlm  
      DO I=1,NS_pot
      if(mod(i,100).eq.0) write(6,*) i
      WRITE(LSTRING,'(i3)') LPOT(I)
      WRITE(MSTRING,'(i3)') MPOT(I)
      OPEN(UNIT=20,FILE='V_lm_DFT_B3LYP/potentiel_l_'//trim(adjustl(lstring)) &
      //'_m_'//trim(adjustl(mstring))//'.dat',STATUS='old')
      DO IR=1,NBR
      READ(20,*) RR(IR),POTlmR,POTlmI
      POTlm(IR,I)=DCMPLX(POTlmR,POTlmI)
      ENDDO
      CLOSE(20)
      ENDDO

! On agrandit la grille en r
      do ir=1,nrini-nbr
      rr(nbr+ir)=rr(nbr)+ir*(REND-rr(nbr))/dfloat(nrini-nbr)
      enddo

! prolongation du potentiel jusqu'a r=RMAX/RK avec le développement multipolaire
      do ir=nbr+1,nrini
      do i=1,NS_pot
      POTlm(ir,i)=rr(nbr)**(lpot(i)+1)*Potlm(nbr,i)/rr(ir)**(lpot(i)+1)
      enddo
      enddo

! couplage entre les états sur la grille en r 
      write(6,*) 'debut du couplage'
      allocate(integrale_3J(NS_pot))
      DO K=1,NS
      if(mod(k,10).eq.0) write(6,*) k
      DO KK=1,NS
      DO I=1,NS_POT
      CALL integrale3ylm(K,KK,i,u)
      integrale_3J(I)=u
      ENDDO
      DO IR=1,NRINI
      V_sum=dcmplx(0.d0,0.d0)
      DO i=1,NS_pot
      IF(MOD(L(k)+L(KK)+LPOT(i),2).EQ.0.AND.&
      abs(L(k)-L(KK)).LE.LPOT(i).AND.&
      L(k)+L(KK).GE.LPOT(i).AND.&
      M(KK)-M(K)+MPOT(i).EQ.0) THEN

      V_sum=V_sum+POTlm(ir,i)*integrale_3J(I)
	endif

      ENDDO
      VV(IR,K,KK)=-2.d0*V_sum
      ENDDO
      ENDDO
      ENDDO


      ALLOCATE(VV_real(NRINI,NS,NS),CVV_real(NRINI,NS,NS),VV_real_bis(NS,NS),VV_bis(NS,NS))
      DO IR=1,NRINI !nbr
      DO K=1,NS
      DO KK=1,NS
      VV_bis(K,KK)=VV(IR,K,KK)
      if(k.eq.kk)then
!      write(60+K,*) rr(ir),VV_bis(k,kk)
      endif
      ENDDO
      ENDDO
      VV_real_bis=MATMUL(C_inv,MATMUL(VV_bis,C))

! verification que Vreel est bien reel (aux erreurs numeriques pres)
      DO K=1,NS
      DO KK=1,NS
       if(dimag(VV_real_bis(K,KK)).GE.(1.D-8)) then
       write(6,*) rr(ir),l(k),m(k),l(kk),m(kk),dreal(VV_real_bis(K,KK)),dimag(VV_real_bis(K,KK))
       stop 259
       endif
      VV_real(IR,K,KK)=VV_real_bis(K,KK)
      ENDDO
      ENDDO
      ENDDO

! preparation du spline
      DO K=1,NS                                                    
      DO KK=1,NS                                                    
      CALL SPLINE(NRINI,RR,VV_real(1,K,KK),CVV_real(1,K,KK),AAA,BBB,CCC)
      ENDDO
      ENDDO

      write(6,*) 'OK'
      CALL NORMNUM

      END

      SUBROUTINE LOGFAC
      USE TABS
      IMPLICIT NONE
      INTEGER :: I
      FAC(1)=0.D0
      DO 10 I=2,300
      FAC(I)=FAC(I-1)+DLOG(DFLOAT(I-1))
  10  CONTINUE
      RETURN
      END
 
      SUBROUTINE integrale3ylm(K,KK,i,u)
!     calcule le couplage entre les états
      USE TABS
      IMPLICIT NONE
      INTEGER :: IR,K,KK,I
      DOUBLE PRECISION :: Z1,Z2
      COMPLEX*16 u

!      CALL LOGFAC

      IF(MOD(L(k)+L(KK)+LPOT(i),2).EQ.0.AND.&
      abs(L(k)-L(KK)).LE.LPOT(i).AND.&
      L(k)+L(KK).GE.LPOT(i).AND.&
      M(KK)-M(K)+MPOT(i).EQ.0) THEN

       CALL CLESPI(2*L(K),2*LPOT(i),2*L(KK),-2*M(K),2*MPOT(i),2*M(KK),Z1)
       CALL CLESPI(2*L(K),2*LPOT(i),2*L(KK),0,0,0,Z2)

      u=z1*z2*(-1.d0)**M(K)*sqrt((2.d0*L(K)+1.d0) &
      *(2.d0*L(KK)+1.d0)*(2.d0*LPOT(i)+1.d0)/(4.d0*PI))
      ELSE
      u=dcmplx(0.d0,0.d0)
      ENDIF
      RETURN
      END

      SUBROUTINE CLESPI(I1,I2,I3,M1,M2,M3,Z)
      USE TABS
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL INIT
      L4=(I1+I2+I3)/2+2
      IF(L4.GT.IDFA) GO TO 1
      IF(M1+M2+M3) 2,3,2        
! ajout par rapport a orion:
    3 IF(M1.EQ.0.AND.M2.EQ.0.AND.M3.EQ.0.AND.MOD(I1+I2+I3,2).NE.0) THEN
      Z=0.D0
      RETURN
      ENDIF
      IZMAX=MIN(I1+I2-I3,I1-M1,I2+M2)/2+1
      IZMIN=MAX(0,I2-I3-M1,I1+M2-I3)/2+1
      IF(IZMAX-IZMIN) 2,4,4       
    4 L1=(I1+I2-I3)/2+1
      L2=(I3+I1-I2)/2+1
      L3=(I3+I2-I1)/2+1
      L5=(I1+M1)/2+1
      L6=(I1-M1)/2+1
      L7=(I2+M2)/2+1
      L8=(I2-M2)/2+1
      L9=(I3+M3)/2+1
      L10=(I3-M3)/2+1
      ABRA=0.5D0*(FAC(L1)+FAC(L2)+FAC(L3)-FAC(L4)+FAC(L5)+FAC(L6)+&
      FAC(L7)+FAC(L8)+FAC(L9)+FAC(L10))   

      K1=(I3-I2+M1)/2+1
      K2=(I3-I1-M2)/2+1
      GROS=250.
      DO 8 II=IZMIN,IZMAX
      I=II-1
      K3=L1-I
      K4=L6-I
      K5=L7-I
      K6=K1+I
      K7=K2+I
      AC(II)=FAC(I+1)+FAC(K3)+FAC(K4)+FAC(K5)+FAC(K6)+FAC(K7)
      IF(AC(II).LT.GROS)  GROS=AC(II)
    8 CONTINUE
      ACCU=0.D0
      SIG=(-1.)**IZMIN
      DO 9 II=IZMIN,IZMAX
      SIG=-SIG
      AC(II)=AC(II)-GROS
      ACCU=ACCU+SIG*DEXP(-AC(II))
    9 CONTINUE
      Z=(-1.)**((I1-I2-M3)/2)*DEXP(ABRA-GROS)*ACCU
! rajout par rapport a orion: zero numerique
      if(dabs(z).le.1.d-9) z=0.d0
      RETURN
   11 WRITE(6,98)
   98 FORMAT(1X,'WARNING: LOGFAC NOT CALLED BEFORE CLESPI')
      STOP
    2 Z=0.D0
      RETURN
    1 WRITE(6,99) 
   99 FORMAT(1X,'LOGSPI: VALUE OF PARAMETER TOO HIGH')
      RETURN
      END


      SUBROUTINE NORMNUM
      USE TABS
!************************************************************************
!  NORMALIZED NUMEROV SELON JOHNSON, J CHEM PHYS 69, 4678 (1978)        *
!************************************************************************
      IMPLICIT NONE
      COMPLEX*16, ALLOCATABLE,DIMENSION(:,:) :: SMAT
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: PSIK
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: Y,V2,U2,R2,PSI2,W2,F2,PSI3,R3,W3,W2G,PSI2_norm
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: FL,GL,FL1,GL1,FLbis,GLbis,GMAT,DELTA
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: PSI_TEST,PSI_TEST_K
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: R
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: COSETA,SINETA,COSETAG,TANETA,SECOND
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FC,GC,FCP,GCP
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: FC1,GC1,FCP1,GCP1
      DOUBLE PRECISION :: AFC,AGC,AFC1,AGC1,BLDENOM,BLNUM,CLNUM,RACT,test1
      INTEGER :: I,J,IR,NN,KOD,lmaxterm,NR_EFF


      ALLOCATE(R(NR))
      ALLOCATE(Y(NS,NS),V2(NS,NS),U2(NS,NS),R2(NS,NS),PSI2(NS,NS),W2(NS,NS),F2(NS,NS))
      ALLOCATE(PSI3(NS,NS),R3(NS,NS),W3(NS,NS),W2G(NS,NS),PSIK(NR,NS,NS),PSI_TEST(NS,NS),PSI_TEST_K(NS,NS),PSI2_norm(NS,NS))
      ALLOCATE(FL(NS,NS),GL(NS,NS),FL1(NS,NS),GL1(NS,NS),FLbis(NS,NS),GLbis(NS,NS),GMAT(NS,NS),DELTA(NS,NS))
      ALLOCATE(COSETA(NS,NS),SINETA(NS,NS),COSETAG(NS,NS),TANETA(NS,NS),SECOND(NS,NS))


      RACT=RSTART

      R2=0.d0

      DO 1 NN=1,NR

!c psi3 sera donc la fonction au pas anterieur; de meme pour les matrices
!c R3 et W3 qui ne sont rien d autre que R2 et W2 au pas anterieur

      PSI3=PSI2
      R3=R2
      W3=W2G

      IF(NN.NE.1) CALL MATINV(NS,R2,DET)    

      RACT=RACT+H
      R(NN)=RACT

      WRITE(6,*) 'NN =',NN,'ract =',ract

      CALL POTM(RACT,V2)
      
!c Calcul de W2 : EQN (23)     
      
      DO I=1,NS
      DO J=1,NS
      W2(I,J)=H*H*V2(I,J)/12.D0
      W2G(I,J)=H*H*V2(I,J)/12.D0
      ENDDO
      W2(I,I)=1.D0+H*H*(RK*RK-L(I)*(L(I)+1.d0)/ract**2.d0+V2(I,I))/12.D0
      W2G(I,I)=1.D0+H*H*(RK*RK-L(I)*(L(I)+1.d0)/ract**2.d0+V2(I,I))/12.D0
      ENDDO

      CALL MATINV(NS,W2,DET)    

!c Calcul de U2 : EQN (24)
      DO I=1,NS
      DO J=1,NS
      U2(I,J)=12.d0*W2(I,J)
      ENDDO
      U2(I,I)=12.d0*W2(I,I)-10.d0
      ENDDO

!c F2 a ete calcule au pas anterieur; je peux terminer le calcul des 
!c fonctions d'onde a ce pas
      IF(NN.GE.2) THEN
      PSI2=MATMUL(W2,F2)
!      IF(mod(NN,5).eq.2) then
      DO I=1,NS
      DO J=1,NS
      PSIk(NN,I,J)=PSI2(I,J)
!      write(6,*) NN,I,J,PSI2(I,J),PSIk(NN,I,J)
!      WRITE(888,*) dsqrt(2.d0/(PI*RK))*PSI2(I,J)/RACT
!      WRITE(2,*) PSI2(I,J)
      ENDDO
      ENDDO
      ENDIF
!      ENDIF

      R2=U2-R2

!c on prepare le calcul des fonctions d onde
      IF(NN.EQ.1) THEN

      DO I=1,NS
      DO J=1,NS
      PSI2(I,J)=0.d0
      ENDDO
      PSI2(I,I)=RACT**(L(I)+1.d0)
      ENDDO

      DO I=1,NS
      DO J=1,NS
      W2(I,J)=H*H*V2(I,J)/12.D0
      W2G(I,J)=H*H*V2(I,J)/12.D0
      ENDDO
      W2(I,I)=1.D0+H*H*(RK*RK-L(I)*(L(I)+1.d0)/ract**2.d0+V2(I,I))/12.D0
      W2G(I,I)=1.D0+H*H*(RK*RK-L(I)*(L(I)+1.d0)/ract**2.d0+V2(I,I))/12.D0
      ENDDO
      F2=MATMUL(W2,PSI2)
      ENDIF

      F2=MATMUL(R2,F2)
! F2 est alors la matrice F au pas suivant

! Calcul de la matrice K (GMAT)
! comme la matrice K est symetrique est que g est diagonale, g.K=K.g
! aussi, les ecritures matricielles f+K.g et f+g.K sont equivalentes
!      go to 1111
      IF(NN.EQ.NR) THEN

      OPEN(UNIT=869,FILE='MATRICE_K',STATUS='new')
! on calcule les fonctions regulieres et irregulieres au dernier pas et 
! a celui d'avant
      ALLOCATE(FC(L(NS)+1),FCP(L(NS)+1),GC(L(NS)+1),GCP(L(NS)+1))
      ALLOCATE(FC1(L(NS)+1),FCP1(L(NS)+1),GC1(L(NS)+1),GCP1(L(NS)+1))

      CALL RCWFN(RK*RACT,-ZEFF/RK,L(1),L(NS),FC,FCP,GC,GCP,1.D-10,KOD)
      CALL RCWFN(RK*(RACT-H),-ZEFF/RK,L(1),L(NS),FC1,FCP1,GC1,GCP1,1.D-10,KOD)

! !!! K est la matrice rangee en vecteurs colonnes pour les Psi's
! !!! i.e. Psi(2,1), par exemple, correspond a f(l(1))*delta(1,2)+K_{21}*g(l(2)), etc....

      DO I=1,NS
      AFC=FC(L(I)+1)
      AGC=GC(L(I)+1)
      AFC1=FC1(L(I)+1)
      AGC1=GC1(L(I)+1)
      DO J=1,NS
      PSI2(I,J)=PSI2(I,J)/RACT
      PSI3(I,J)=PSI3(I,J)/(RACT-H)
      BLNUM=PSI2(I,J)*agc1/(ract-H)-PSI3(I,J)*agc/ract
      CLNUM=PSI2(I,J)*afc1/(ract-H)-PSI3(I,J)*afc/ract
      BLDENOM=-AFC*AGC1+AFC1*AGC
      BLDENOM=BLDENOM/((RACT-H)*RACT)
      SINETA(I,J)=CLNUM/BLDENOM
      COSETA(I,J)=-BLNUM/BLDENOM
      COSETAG(I,J)=COSETA(I,J)
      FL1(I,J)=0.d0
      GL1(I,J)=0.d0
      FL(I,J)=0.d0
      GL(I,J)=0.d0
      ENDDO
      FL(I,I)=FC(L(I)+1)
      GL(I,I)=GC(L(I)+1)
      FL1(I,I)=FC1(L(I)+1)
      GL1(I,I)=GC1(L(I)+1)
      ENDDO

      CALL MATINV(NS,COSETA,DET)

!c  on calcule la matrice K a partir des matrices des fonctions regulieres 
!c  et irregulieres : EQN (A16) (A17) et (A19)
      FL1=MATMUL(W3,FL1)
      GL1=MATMUL(W3,GL1)
      FL=MATMUL(W2G,FL)
      GL=MATMUL(W2G,GL)
      FL1=MATMUL(R3,FL1)
      GL1=MATMUL(R3,GL1)

      FL1=FL1-FL
      GL1=GL1-GL

      CALL MATINV(NS,GL1,DET)
      GMAT=MATMUL(GL1,FL1)

      GMAT=-GMAT ! le (-) vient de numerov

      DEALLOCATE(FC,FCP,GC,GCP)
      DEALLOCATE(FC1,FCP1,GC1,GCP1)

!c comparaison entre la matrice K et B.A-1 (A=coseta B=sineta)
!c alors que GMAT est K, COSETA est B.A-1

      TANETA=MATMUL(SINETA,COSETA)

      WRITE(6,*) 'KMATRIX'
      DO I=1,NS
      DO J=1,NS
	write(66,*) I,J,SINETA(I,J),COSETA(I,J),COSETAG(I,J)
      WRITE(6,52) L(J),M(J),L(I),M(I),TANETA(J,I),GMAT(J,I)
 52   FORMAT(' COMP: L"=',I3,' M"=',I3,' L=',I3,' M=',I3,' B.A-1=',D13.6,' GMAT=',D13.6)
      WRITE(869,*) GMAT(I,J)
      ENDDO
      ENDDO
      write(6,*) 'ok'
!
      DO I=1,NS
      WRITE(6,53) L(I),M(I),GMAT(I,I),ATAN(GMAT(I,I))/pi
 53   FORMAT(I6,I6,' GMAT=',D13.6,' DEPHASAGE=',D13.6)
      ENDDO
!
      CLOSE(869)
      ENDIF

 1    CONTINUE


      NR_EFF=0
      OPEN(UNIT=888,FILE='PSI_K',STATUS='new')
      DO IR=2,NR,5 
      NR_EFF=NR_EFF+1
      ENDDO
      WRITE(888,*) NR_EFF

      DO IR=2,NR,5 
      WRITE(888,*) R(IR)
      RACT=R(IR)
      DO I=1,NS
      DO J=1,NS
      psi2(I,J)=PSIk(ir,I,J)
      ENDDO
      ENDDO
      PSI2_norm=MATMUL(PSI2,COSETA)                                         !COSETA est en fait COSETA-1
      DO I=1,NS
      DO J=1,NS
      WRITE(888,*) dsqrt(2.d0/(PI*RK))*PSI2_norm(I,J)/RACT
      ENDDO
      ENDDO
      ENDDO
      CLOSE(888)

      RETURN
      END

      SUBROUTINE POTM(RACT,V)
!************************************************************************
!* CONSTRUCT DIABATIC POTENTIAL MATRIX AT R=RACT                        *
!* SPLINE INTERPOLATION FOR ENERGIES AND COUPLING                       *
!************************************************************************
      USE TABS
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(NS,NS) :: V
      DOUBLE PRECISION :: RACT,SPL,R_internuc,Z1,Z2,THETA,PHI,PLGNDR
      INTEGER :: I,J,ipot,jpot,m3,IDF
      DOUBLE PRECISION,DIMENSION(101) :: facto_simple
      COMPLEX*16,DIMENSION(10000) :: YLM
      COMPLEX*16 POT_lm,u
      COMPLEX*16,DIMENSION(NS,NS) :: VCOMPLEX

      DO I=1,NS
      DO J=1,I
        V(I,J)=SPL(NRINI,RR,VV_real(1,I,J),CVV_real(1,I,J),RACT)
        V(J,I)=V(I,J)

!      R0=DSQRT(L(I)*(L(I)+1.D0))/RK
!      IF(RACT.LE.R0+(LOG(1.D-2)/RK).AND.I.NE.J) THEN
!      V(I,J)=0.D0
!      V(J,I)=0.D0
!      ENDIF

!	IF(I.eq.J)THEN
!        V(I,J)=2.d0/RACT
!        V(J,I)=V(I,J)
!	ELSE
!	V(I,J)=0.d0
!        V(J,I)=V(I,J)
!	ENDIF

!       V(I,J)=0.D0
!       V(J,I)=0.D0

      ENDDO
      ENDDO

      RETURN                                                           
      END                                                              


      SUBROUTINE SPLINE(N,X,Y,CM,ALPHA,BETA,B)                        
      IMPLICIT REAL*8(A-H,O-Z)                                         
      DIMENSION ALPHA(N),BETA(N),B(N)                                
      DIMENSION X(N),Y(N),CM(N)                                       
      COMMON/KODSPL/KOD,K,IW,ILIGNE
      IW=6                                                            
      ILIGNE=5                                                        
      K=2                                                            
      CM(1)=0.                                                        
      CM(N)=0.                                                         
      KOD=0                                                          
      C=X(2)-X(1)                                                    
      E=(Y(2)-Y(1))/C                                               
      DO 10 I=3,N                                                   
      I2=I-2                                                        
      A=C                                                             
      C=X(I)-X(I-1)                                                   
      IF(A*C.GT.0.0) GO TO 1                                          
      KOD=-1                                                          
      IF(IW.GT.0) WRITE(IW,91) X(I-1),X(I)                            
   91 FORMAT(1X,'ERROR IN THE DATA FOR SPLINE'/' CHECK ABCISSA',E20.8,10 &
      X,'AND',E20.8,/' SPLINE STOPPED')                              
      RETURN                                                         
    1 ALPHA(I2)=(A+C)/3.                                              
      BETA(I2)=C/6.                                                
      D=E                                                             
      E=(Y(I)-Y(I-1))/C                                           
   10 B(I2)=E-D                                                      
      CALL TRIDIA(ALPHA,BETA,BETA,B,CM(2),N-2)                       
      K=2                                                           
      RETURN                                                         
      END                                                          
!C                                                                  
      SUBROUTINE TRIDIA(ALPHA,BETA,GAMMA,B,X,N)                     
      IMPLICIT REAL*8(A-H,O-Z)                                         
      DIMENSION ALPHA(N),BETA(N),GAMMA(N),B(N),X(N)                  
      DO 10 I=2,N                                                    
      RAP=BETA(I-1)/ALPHA(I-1)                                        
      ALPHA(I)=ALPHA(I)-RAP*GAMMA(I-1)                               
   10 B(I)=B(I)-RAP*B(I-1)                                              
      X(N)=B(N)/ALPHA(N)                                               
      DO 20 J=2,N                                                     
      I=N-J+1                                                         
   20 X(I)=(B(I)-GAMMA(I)*X(I+1))/ALPHA(I)                            
      RETURN                                                          
      END                                                             
!C                                                                    
      FUNCTION SPL(N,X,Y,M,T)                                         
      IMPLICIT REAL*8(A-H,O-Z)                                        
      DIMENSION X(N),Y(N)                                              
      REAL*8 M(N)                                                      
      COMMON/KODSPL/KOD,K,IW,ILIGNE
    1 FORMAT(1X,'SPL: TENTATIVE OF EXTRAPOLATION- ABCISSA=',1PE18.8)  
      IF(X(2).LT.X(1)) GO TO 100                                      
      IF(T.LT.X(1)) GO TO 30                                          
      IF(T.GT.X(N)) GO TO 40                                          
      IF(K.GT.N) K=2                                                  
      IF(T.LE.X(K)) GO TO 11                                          
   10 K=K+1                                                           
      IF(T.GT.X(K)) GO TO 10                                          
      GO TO 20                                                        
   11 IF(T.GE.X(K-1)) GO TO 20                                        
      K=K-1                                                           
      GO TO 11                                                        
  100 IF(T.GT.X(1)) GO TO 30                                          
      IF(T.LT.X(N)) GO TO 40                                          
      IF(K.GT.N) K=2                                                  
      IF(T.GE.X(K)) GO TO 111                                           
  110 K=K+1                                                             
      IF(T.LT.X(K)) GO TO 110                                           
      GO TO 20                                                          
  111 IF(T.LE.X(K-1)) GO TO 20                                          
      K=K-1                                                             
      GO TO 111                                                         
   20 CONTINUE                                                          
      F=X(K)-T                                                          
      G=T-X(K-1)                                                        
      E=F+G                                                             
      SPL=(-G*F*(M(K-1)*(F+E)+M(K)*(G+E))+6.*(G*Y(K)+F*Y(K-1)))/(6.*E)  
      KOD=0                                                             
      RETURN                                                            
   30 E=X(2)-X(1)                                                       
      SPL=((Y(2)-Y(1))/E-M(2)*E/6.)*(T-X(1))+Y(1)                       
      KOD=1                                                             
      WRITE(6,1) T  
      RETURN                                                           
   40 E=X(N)-X(N-1)                                                    
      SPL=((Y(N)-Y(N-1))/E+M(N-1)*E/6.)*(T-X(N))+Y(N)                  
      KOD=2                                                            
      WRITE(6,1) T                     
      RETURN                                                            
      END                                                               
!c
      SUBROUTINE MATINV(NORDER,A,DET)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NORDER,NORDER),IK(NORDER),JK(NORDER) 
!C INVERT A SYMMETRIC MATRIX AND CALCULATE ITS DETERMINANT.             
!C        DESCRIPTION OF PARAMETERS:                                    
!C A  - INPUT MATRIX WHICH IS REPLACED BY ITS INVERSE                   
!C DET- DETERMINANT OF MATRIX                                           
      DET=1.0                                                           
      DO 100 K=1,NORDER                                                 
!C FIND THE LARGEST ELEMENT A(I,J) IN REST MATRIX                       
      AMAX=0.0                                                          
 21   DO 30 I=K,NORDER                                                  
      DO 29 J=K,NORDER                                                  
      IF(DABS(AMAX).GT.DABS(A(I,J)))GOTO 29                             
      AMAX=A(I,J)                                                       
      IK(K)=I                                                           
      JK(K)=J                                                          
 29   CONTINUE                                                          
 30   CONTINUE                                                        
!C INTERCHANGE ROWS AND COLUMS TO PUT AMAX IN A(K,K)                 
      IF(AMAX.NE.0.0)GOTO 41                                        
      DET=0.0                                                      
      GOTO 140                                                    
 41   I=IK(K)                                                    
      IF(I-K)21,51,43                                           
 43   DO 50 J=1,NORDER                                         
      SAVE=A(K,J)                                           
      A(K,J)=A(I,J)                                           
 50   A(I,J)=-SAVE                                           
 51   J=JK(K)                                              
      IF(J-K)21,61,53                                     
 53   DO 60 I=1,NORDER                                   
      SAVE=A(I,K)                                                      
      A(I,K)=A(I,J)                                                   
 60   A(I,J)=-SAVE                                                   
!C ACCUMULATE ELEMENTS OF INVERSE MATRIX                            
 61   DO 70 I=1,NORDER                                             
      IF(I.EQ.K)GOTO 70                                           
      A(I,K)=-A(I,K)/AMAX                                        
 70   CONTINUE                                                  
      DO 80 I=1,NORDER                                         
      DO 80 J=1,NORDER                                        
      IF(I.EQ.K)GOTO 80                                      
      IF(J.EQ.K)GOTO 80                                     
      A(I,J)=A(I,J)+A(I,K)*A(K,J)                          
 80   CONTINUE                                            
      DO 90 J=1,NORDER                                   
      IF(J.EQ.K)GOTO 90                                 
      A(K,J)=A(K,J)/AMAX                             
 90   CONTINUE                                         
      A(K,K)=1.0/AMAX                                 
 100  DET=DET*AMAX                                  
!C RESTORE ORDERING OF MATRIX                      
      DO 130 L=1,NORDER                           
      K=NORDER-L+1                               
      J=IK(K)                                   
      IF(J.LE.K)GOTO 111                       
      DO 110 I=1,NORDER                       
      SAVE=A(I,K)                            
      A(I,K)=-A(I,J)                        
 110  A(I,J)=SAVE                          
 111  I=JK(K)                             
      IF(I.LE.K)GOTO 130                 
      DO 120 J=1,NORDER                 
      SAVE=A(K,J)                     
      A(K,J)=-A(I,J)                   
 120  A(I,J)=SAVE                 
 130  CONTINUE                       
 140  RETURN                       
      END                           
                                 
      SUBROUTINE RCWFN(RHO,ETA,MINL0,MAXL,FC,FCP,GC,GCP,ACCUR,KOD)      
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      REAL*8 K                                                          
      DIMENSION FC(MAXL+1),FCP(MAXL+1),GC(MAXL+1),GCP(MAXL+1)
      DATA OVRFLW/1.0D+60/,UNDFLW/1.D-36/                               
!      FLOAT(I)=DFLOAT(I)                                                
!      EXP(X)=DEXP(X)                                                    
!      ALOG(X)=DLOG(X)                                                   
!      ABS(X)=DABS(X)                                                    
!      IFIX(X)=IDINT(X)                                                  
!      SQRT(X)=DSQRT(X)                                                  
! *** COULOMB WAVEFUNCTIONS CALCULATED AT R = RHO BY THE                
! *** CONTINUED-FRACTION METHOD OF STEED   MINL,MAXL ARE ACTUAL L-VALUES
! *** SEE BARNETT FENG STEED AND GOLDFARB COMPUTER PHYSICS COMMUN 1974  
      KOD=0                                                             
      MINL=IABS(MINL0)                                                  
      IBORNE=0                                                          
      IF(MINL0.LT.0) IBORNE=MINL                                        
      ACC  = ACCUR                                                      
      IF(ACC.GT.1.0E-6) ACC=1.0E-6                                      
      IF(ACC.LT.1.0E-14) ACC=1.0E-14                                    
      R    = RHO                                                        
      KTR  = 1                                                          
      LMAX = MAXL                                                       
      LMIN1= MINL + 1                                                   
      XLL1 = FLOAT(MINL*LMIN1)                                          
      ETA2 = ETA*ETA                                                    
      TURN = ETA + SQRT(ETA2 + XLL1)                                    
      IF(R.LT.TURN) KTR=-1                                              
      KTRP = KTR                                                        
      GO TO 2                                                           
1     R    = TURN                                                       
      TF   = F                                                          
      TFP  = FP                                                         
      LMAX = MINL                                                       
      KTRP = 1                                                          
2     ETAR = ETA*R                                                      
      RHO2 =   R*R                                                      
      PL   = FLOAT(LMAX + 1)                                            
      PMX  = PL + 0.5                                                   
! *** CONTINUED FRACTION FOR FP(MAXL)/F(MAXL)  XL IS F  XLPRIME IS FP **
      FP  = ETA/PL + PL/R                                               
      DK  = ETAR*2.0                                                    
      DEL = 0.0                                                         
      D   = 0.0                                                         
      F   = 1.0                                                         
      K   = (PL*PL - PL + ETAR)*(2.0*PL - 1.0)                          
      IF(DABS(PL*PL+PL+ETAR).GT.UNDFLW) GO TO 3                         
1000  R   = R + ACC                                                     
      GO TO 2                                                           
3     H   = (PL*PL + ETA2)*(1.0 - PL*PL)*RHO2                           
      K   = K + DK + PL*PL*6.0                                          
      D   =  1.0/(D*H + K)                                              
      DEL =  DEL*(D*K - 1.0)                                            
      IF(PL.LT.PMX) DEL = -R*(PL*PL + ETA2)*(PL + 1.0)*D/PL             
      PL  = PL + 1.0                                                    
      FP  = FP + DEL                                                    
      IF(D.LT.0.0) F = -F                                               
      IF(PL.GT.20000.) GOTO 102                                         
      IF(DABS(FP).LT.UNDFLW) GO TO 1000                                 
      IF(ABS(DEL/FP).GE.ACC) GO TO 3                                    
      FP  = F*FP                                                        
      IF( LMAX.EQ.MINL) GO TO 5                                         
      IPIC=LMAX+1-IBORNE                                                
      FC (IPIC) = F                                                     
      FCP(IPIC) = FP                                                    
! *** DOWNWARD RECURSION TO MINL FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      L  = LMAX                                                         
      DO 4 LP  = LMIN1,LMAX                                             
      PL = FLOAT(L)                                                     
      IPIC=L-IBORNE                                                     
      GC (IPIC+1) = ETA/PL + PL/R                                       
      GCP(IPIC+1) = SQRT(ETA2 + PL*PL)/PL                               
      FC (IPIC)   = (GC(IPIC+1)*FC(IPIC+1) + FCP(IPIC+1))/GCP(IPIC+1)   
      FCP(IPIC)   =  GC(IPIC+1)*FC(IPIC)   - GCP(IPIC+1)*FC(IPIC+1)     
      IF(ABS(FCP(IPIC)).GT.OVRFLW) GOTO 103                             
4     L  = L - 1                                                        
      IPIC=LMIN1-IBORNE                                                 
      F  = FC (IPIC)                                                    
      FP = FCP(IPIC)                                                    
5     IF(KTRP.EQ.-1) GO TO 1                                            
! *** REPEAT FOR R = TURN IF RHO LT TURN                                
! *** NOW OBTAIN P + I.Q FOR MINL FROM CONTINUED FRACTION (32)          
! *** REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING REAL*8      
      P  = 0.0                                                          
      Q  = R - ETA                                                      
      PL = 0.0                                                          
      AR = -(ETA2 + XLL1)                                               
      AI =   ETA                                                        
      BR = 2.0*Q                                                        
      BI = 2.0                                                          
      WI = 2.0*ETA                                                      
      DR =   BR/(BR*BR + BI*BI)                                         
      DI =  -BI/(BR*BR + BI*BI)                                         
      DP = -(AR*DI + AI*DR)                                             
      DQ =  (AR*DR - AI*DI)                                             
6     P  =  P + DP                                                      
      Q  =  Q + DQ                                                      
      PL = PL + 2.0                                                     
      AR = AR + PL                                                      
      AI = AI + WI                                                      
      BI = BI + 2.0                                                     
      D  = AR*DR - AI*DI + BR                                           
      DI = AI*DR + AR*DI + BI                                           
      T  = 1.0/(D*D + DI*DI)                                            
      DR =  T*D                                                         
      DI = -T*DI                                                        
      H  = BR*DR - BI*DI - 1.0                                          
      K  = BI*DR + BR*DI                                                
      T  = DP*H  - DQ*K                                                 
      DQ = DP*K  + DQ*H                                                 
      DP = T                                                            
      IF(PL.GT.46000.) GOTO 101                                         
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6                
      P  = P/R                                                          
      Q  = Q/R                                                          
! *** SOLVE FOR FP,G,GP AND NORMALISE F  AT L=MINL                      
      G  = (FP - P*F)/Q                                                 
      GP = P*G - Q*F                                                    
      FPF=FP/F                                                          
      GF =G/F                                                           
      GPF=GP/F                                                          
      W=1.0/(ABS(F)*SQRT(FPF*GF-GPF))                                   
      G  = W*G                                                          
      GP = W*GP                                                         
      IF(KTR.EQ.1) GO TO 8                                              
      F  = TF                                                           
      FP = TFP                                                          
      LMAX = MAXL                                                       
! *** INTEGRATION OF G(MINL) AND GP(MINL) A R DECROISSANT DEPUIS TURN   
! *** PAR UNE METHODE D'EXTRAPOLATION A PAS VARIABLE.                   
! *** CF. DIPLOMARBEIT VON HANS-GUNTER HUSSELS,                         
! ***     AIMABLEMENT COMMUNIQUE PAR LE PROF. DR. R BURLIRSCH.          
      X=TURN                                                            
      Y=G                                                               
      DYDX=GP                                                           
      H=.1                                                              
      NSTEP=0                                                           
   21 DELTA=RHO-X                                                       
      DELABS=DABS(DELTA)                                                
      IF(DELABS/RHO.LT.1.D-14) GOTO 22                                  
      H=DSIGN(DMIN1(DABS(H),DELABS),DELTA)                              
      CALL DIFSYL(ETA+ETA,XLL1,ACC,H,X,Y,DYDX)                          
      IF(ABS(DYDX).GT.OVRFLW) GOTO 105                                  
      IF(H.EQ.0.0) GOTO 107                                             
      NSTEP=NSTEP+1                                                     
      IF(NSTEP.GT.20000) GOTO 106                                       
      GOTO 21                                                           
   22 G=Y                                                               
      GP=DYDX                                                           
      IF(SQRT(ABS(FP))*SQRT(ABS(GP)).GT.SQRT(OVRFLW)) GOTO 104          
      W  = 1.0/(FP*G - F*GP)                                            
! *** UPWARD RECURSION FROM GC(MINL) AND GCP(MINL),STORED VALUES ARE R,S
! *** RENORMALISE FC,FCP FOR EACH L-VALUE                               
8     IPIC=LMIN1-IBORNE                                                 
      GC (IPIC) = G                                                     
      GCP(IPIC) = GP                                                    
      IF(LMAX.EQ.MINL) GO TO 10                                         
      DO  9  L = LMIN1,LMAX                                             
      IPIC=L-IBORNE                                                     
      T        = GC(IPIC+1)                                             
      GC (IPIC+1) = (GC(IPIC)*GC (IPIC+1) - GCP(IPIC))/GCP(IPIC+1)      
      GCP(IPIC+1) =  GC(IPIC)*GCP(IPIC+1) - GC(IPIC+1)*T                
      IF(ABS(GCP(IPIC+1)).GT.OVRFLW) GOTO 103                           
      FC (IPIC+1) = W*FC (IPIC+1)                                       
9     FCP(IPIC+1) = W*FCP(IPIC+1)                                       
      IPIC=LMIN1-IBORNE                                                 
      FC (IPIC) = FC (IPIC)*W                                           
      FCP(IPIC) = FCP(IPIC)*W                                           
      RETURN                                                            
10    IPIC=LMIN1-IBORNE                                                 
      FC (IPIC) = W*F                                                   
      FCP(IPIC) = W*FP                                                  
      RETURN                                                            
11    W  = 0.0                                                          
      G  = 0.0                                                          
      GP = 0.0                                                          
      GO TO 8                                                           
  101 KOD=1                                                             
      GOTO 11                                                           
  102 KOD=2                                                             
      GOTO 11                                                           
  103 KOD=3                                                             
      GOTO 11                                                           
  104 KOD=4                                                             
      GOTO 11                                                           
  105 KOD=5                                                             
      GOTO 11                                                           
  106 KOD=6                                                             
      GOTO 11                                                           
  107 KOD=7                                                             
      GOTO 11                                                           
      END      
                                        
      SUBROUTINE DIFSYL(ETA2,XLL1,EPS,H,X,Y,DYDX)                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      REAL*4 FA,FV,ETA,FY,EP,FS,FLOAT,SNGL                              
      DIMENSION YR(2),YS(2),DT(2,7),D(7),S(2),EP(4)                     
      LOGICAL   KONV,BO,KL,GR                                           
      DATA EP/0.4E-1,0.16E-2,0.64E-4,0.256E-5/                          
      DATA N/1/                                                         
      COU(X)=ETA2/X+XLL1/(X*X)-1.D0                                     
      JTI=0                                                             
      FY=1.                                                             
      ETA=SNGL(DABS(EPS))                                               
      IF(ETA.LT.1.E-11) ETA=1.E-11                                      
      DZ2=COU(X)*Y                                                      
   10 XN=X+H                                                            
      BO=.FALSE.                                                        
      S(1)=DABS(Y)                                                      
      S(2)=DABS(DYDX)                                                   
      M=1                                                               
      JR=2                                                              
      JS=3                                                              
      DO 260 J=1,10                                                     
      IF(.NOT.BO) GOTO 200                                              
      D(2)=1.777777777777778D0                                          
      D(4)=7.111111111111111D0                                          
      D(6)=2.844444444444444D1                                          
      GOTO 201                                                          
  200 D(2)=2.25D0                                                       
      D(4)=9.D0                                                         
      D(6)=3.6D1                                                        
  201 IF(J.LE.7) GOTO 202                                               
      L=7                                                               
      D(7)=6.4D1                                                        
      GOTO 203                                                          
  202 L=J                                                               
      D(L)=M*M                                                          
  203 KONV=L.GT.3                                                       
      B=H/DFLOAT(M)                                                     
      G=B*0.5D00                                                        
      YS(2)=DYDX+G*DZ2                                                  
      YS(1)=Y+B*YS(2)                                                   
      M=M-1                                                             
      IF(IABS(M).LT.1) GOTO 221                                         
      DO 220 K=1,M                                                      
      DY2=COU(X+DFLOAT(K)*B)*YS(1)                                      
      I=1                                                               
      U=DABS(YS(I))                                                     
      IF(U.GT.S(I)) S(I)=U                                              
      NI=2                                                              
      V=YS(NI)+B*DY2                                                    
      YS(I)=YS(I)+B*V                                                   
      YS(NI)=V                                                          
      V=DABS(V)                                                         
      IF(V.GT.S(NI)) S(NI)=V                                            
  220 CONTINUE                                                          
  221 DY2=COU(XN)*YS(1)                                                 
      YS(2)=YS(2)+G*DY2                                                 
      KL=L.LT.2                                                         
      GR=L.GT.5                                                         
      FS=0.                                                             
      M=N+N                                                             
      DO 233 I=1,M                                                      
      V=DT(I,1)                                                         
      C=YS(I)                                                           
      U=DABS(C)                                                         
      IF(U.GT.S(I)) S(I)=U                                              
      DT(I,1)=C                                                         
      TA=C                                                              
      IF(KL) GOTO 233                                                   
      DO 231 K=2,L                                                      
      B1=D(K)*V                                                         
      B=B1-C                                                            
      W=C-V                                                             
      U=V                                                               
      IF(B.EQ.0.D0) GOTO 230                                            
      B=W/B                                                             
      U=C*B                                                             
      C=B1*B                                                            
  230 V=DT(I,K)                                                         
      DT(I,K)=U                                                         
  231 TA=U+TA                                                           
      IF(.NOT.KONV) GOTO 232                                            
      IF(DABS(YR(I)-TA).GT.S(I)*ETA) KONV=.FALSE.                       
  232 IF(GR.OR.S(I).EQ.0.D0) GOTO 233                                   
      FV=DABS(W)/S(I)                                                   
      IF(FS.LT.FV) FS=FV                                                
  233 YR(I)=TA                                                          
      IF(FS.EQ.0.D0) GOTO 250                                           
      FA=FY                                                             
      K=L-1                                                             
      FY=(EP(K)/FS)**(1./FLOAT(L+K))                                    
      IF(L.EQ.2) GOTO 240                                               
      IF(FY.LT.0.7*FA) GOTO 250                                         
  240 IF(FY.GT.0.7) GOTO 250                                            
      H=H*FY                                                            
      JTI=JTI+1                                                         
      IF(JTI.GT.5) GOTO 30                                              
      GOTO 10                                                           
  250 IF(KONV) GOTO 20                                                  
      D(3)=4.D0                                                         
      D(5)=1.6D1                                                        
      BO=.NOT.BO                                                        
      M=JR                                                              
      JR=JS                                                             
  260 JS=M+M                                                            
      H=H*0.5D0                                                         
      GOTO 10                                                           
   20 X=XN                                                              
      H=H*FY                                                            
      Y=YR(1)                                                           
      DYDX=YR(2)                                                        
      RETURN                                                            
   30 H=0.D00                                                           
      RETURN                                                            
      END          
!c
      SUBROUTINE CMATINV(NORDER,A,DETC)                            
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A,DETC,AMAX,SAVE
      DIMENSION A(NORDER,NORDER),IK(NORDER),JK(NORDER)
!C INVERT A GENERAL COMPLEX MATRIX AND CALCULATE ITS DETERMINANT.       
!C AT PRESENT NDIM MUST BE LESS THAN TEN                               
      DETC=CMPLX(1.0,0.0)                                             
      DO 100 K=1,NORDER                                              
!C FIND THE LARGEST ELEMENT A(I,J) IN REST MATRIX                   
      AMAX=CMPLX(0.0,0.0)                                          
 21   DO 30 I=K,NORDER                                            
      DO 29 J=K,NORDER                                           
      IF(CDABS(AMAX).GT.CDABS(A(I,J))) GOTO 29                  
      AMAX=A(I,J)                                              
      IK(K)=I                                                 
      JK(K)=J                                                
 29   CONTINUE                                              
 30   CONTINUE                                             
!C INTERCHANGE ROWS AND COLUMS TO PUT AMAX IN A(K,K)      
      IF(AMAX.NE.DCMPLX(0.D0,0.D0))GOTO 41               
      DETC=0.0                                          
      GOTO 140                                         
 41   I=IK(K)                                         
      IF(I-K)21,51,43                                
 43   DO 50 J=1,NORDER                              
      SAVE=A(K,J)                                  
      A(K,J)=(A(I,J))                             
 50   A(I,J)=-(SAVE)                             
 51   J=JK(K)                                   
      IF(J-K)21,61,53                          
 53   DO 60 I=1,NORDER                        
      SAVE=A(I,K)                            
      A(I,K)=(A(I,J))                       
 60   A(I,J)=-(SAVE)                       
!C ACCUMULATE ELEMENTS OF INVERSE MATRIX  
 61   DO 70 I=1,NORDER                   
      IF(I.EQ.K)GOTO 70                 
      A(I,K)=-A(I,K)/AMAX              
 70   CONTINUE                        
      DO 80 I=1,NORDER               
      DO 80 J=1,NORDER              
      IF(I.EQ.K)GOTO 80            
      IF(J.EQ.K)GOTO 80           
      A(I,J)=A(I,J)+A(I,K)*A(K,J)
 80   CONTINUE                  
      DO 90 J=1,NORDER         
      IF(J.EQ.K)GOTO 90       
      A(K,J)=A(K,J)/AMAX     
 90   CONTINUE              
      A(K,K)=1.0/AMAX      
 100  DETC=DETC*AMAX
!C RESTORE ORDERING OF MATRIX                                          
      DO 130 L=1,NORDER                                               
      K=NORDER-L+1                                                   
      J=IK(K)                                                       
      IF(J.LE.K)GOTO 111                                           
      DO 110 I=1,NORDER                                           
      SAVE=A(I,K)                                                
      A(I,K)=-(A(I,J))                                          
 110  A(I,J)=(SAVE)                                           
 111  I=JK(K)                                                
      IF(I.LE.K)GOTO 130                                    
      DO 120 J=1,NORDER                                    
      SAVE=A(K,J)                                         
      A(K,J)=-(A(I,J))                                   
 120  A(I,J)=(SAVE)                                     
 130  CONTINUE                                         
 140  RETURN                                          
      END                                            

        double precision function plgndr(l,m,x)
        implicit none
        integer l,m,i,ll
        double precision x,fact,pll,pmm,pmmp1,somx2
        if(m.LT.0.OR.m.GT.l.OR.abs(x).GT.1) stop 'bad arguments in plgndr'
        pmm=1.d0
        if (m.gt.0) then
                somx2=sqrt((1.d0-x)*(1.d0+x))
                fact=1.d0
                do i=1,m
                        pmm=-pmm*fact*somx2
                        fact=fact+2.d0
                enddo
        endif
        if(l.EQ.m) then
                plgndr=pmm
        else
                pmmp1=x*(2*m+1)*pmm
                if(l.EQ.m+1) then
                        plgndr=pmmp1
                else
                        do ll=m+2,l
                                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                                pmm=pmmp1
                                pmmp1=pll
                        enddo
                        plgndr=pll
                endif
        endif
        return
        end function plgndr

