!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! calculo da area de um circulo de raio R = 3

! A = (\int _{0}^{2 \pi} d\phi) (\int_{0}^{R}) r dr = \pi R^{2} 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  


	implicit double precision (a-z) 

	integer ndim, NPRN, IGRAPH, NCALL1, ITMX1      	
	common/limites/ ph_i, ph_f, r_i, r_f
        COMMON/FER/AVGI     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        COMMON/cisma/NDIM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! necessarios para a VEGAS
        external AINTEG     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	common/raio/ Rc
	external phint !integrando chamado como funcao external


	pi = dacos(-1.d0)
	
	Rc = 3.d0 ! escolha do raio do circulo



C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!	LIMITES DE INTEGRACAO          
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	ph_i = 0.d0
	ph_f = 2.d0*pi 
	
	r_i = 0.d0
	r_f = Rc


C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  
!	CHAMA A VEGAS          
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

      NDIM = 2			! numero de dimensoes
      ACC = 1.D-16              ! precisao
      NPRN = 0 			! 0 ou 1
      IGRAPH = 0
      NCALL1 = 2000 ! 1000	! numero de chamadas
      ITMX1 = 20 ! 20		! numero de iteracoes
C     
      AVGI = 0.D0
      CALL VEGAS(AINTEG,ACC,NDIM,NCALL1,ITMX1,NPRN,IGRAPH)
      RESU = AVGI         
	

C     ------------------------------------------------------
C     REFINAMENTO DOS RESULTADOS !!!!!
      NCALL1 = 20000!10000
      ITMX1 = 30! 30
C     
      AVGI = 0.D0
      CALL VEGAS1(AINTEG,ACC,NDIM,NCALL1,ITMX1,NPRN,IGRAPH)


      RESU = AVGI ! RESULTADO          

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	INTEGRAL = RESU




	write(*,*) 'resultado numerico: ',integral
	write(*,*) 'area analitico: ', pi * Rc**2.d0

	end





!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Integracao com a VEGAS

	DOUBLE PRECISION function AINTEG(ZZ)  ! funcao onde definimos o integrando
        IMPLICIT DOUBLE PRECISION(A-Z)
       	common/limites/ ph_i, ph_f, r_i, r_f
       	integer ndim, NPRN, IGRAPH, NCALL1, ITMX1      	
        COMMON/FER/AVGI
        COMMON/cisma/NDIM
 

      DIMENSION ZZ(NDIM)! define o tamanho do vetor
C     
C
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C DEFINICOES DE VARIAVEIS EM TERMOS DOS LIMITES DE INTEGRACAO
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	ph = ph_i + (ph_f - ph_i) * ZZ(1)      
			
                                 
	r = r_i + (r_f - r_i) * ZZ(2)      
                        


C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C	JACOBIANO da rotina
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	dph = ph_f - ph_i

	dr = r_f - r_i

	JAC = dph * dr
	
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



	Test = r


C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C	INTEGRANDO
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



	if(Test .lt. 1.d-150 ) then


		AINTEG = 0.d0
!	
			! JAC = chamado de jacobiano da VEGAS. Igual ao produto da variação dos limites de integracao de todas as 										variaveis
	else

			AINTEG = Test * JAC

		!       integrando * jacobiano da VEGAS (nao eh a mesma coisa que o jacobiano da integral)  

	end if
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      RETURN
      END	














CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Rotina de Integracao
c
c     sgs0(li,ls,precisao,integrando)
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION SGS0 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS8,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS8 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS8 ( F,U,C )
      SG = SGS8 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS0=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS8 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS8 = S * H
      RETURN
      END
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS2 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS6,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS6 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS6 ( F,U,C )
      SG = SGS6 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS2=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS6 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS6 = S * H
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS4 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS9,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS9 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS9 ( F,U,C )
      SG = SGS9 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS4=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS9 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS9 = S * H
      RETURN
      END
c
c



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS5 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS10,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS10 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS10 ( F,U,C )
      SG = SGS10 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS5=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS10 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS10 = S * H
      RETURN
      END
c
c



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS3 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS7,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS7 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS7 ( F,U,C )
      SG = SGS7 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS3=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS7 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS7 = S * H
      RETURN
      END
c
c


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS11 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS12,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS12 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS12 ( F,U,C )
      SG = SGS12 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS11=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS12 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS12 = S * H
      RETURN
      END
c
c



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS21 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS22,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS22 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS22 ( F,U,C )
      SG = SGS22 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS21=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS22 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS22 = S * H
      RETURN
      END
c
c


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      DOUBLE PRECISION FUNCTION SGS23 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS24,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS24 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS24 ( F,U,C )
      SG = SGS24 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS23=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS24 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS24 = S * H
      RETURN
      END




C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C


C    VEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGASVEGAS
C   20/11/89 911201938  MEMBER NAME  HBVEGAS  (AXO)         FVS
      SUBROUTINE VGPDAT(A)                                              00000100
C  Dummy routine for VGPDAT.                                            00000110
C                                                                       00000120
      DUMMY=A                                                           00000130
      RETURN                                                            00000200
      END                                                               00000300
      SUBROUTINE VGDAT
C
C  INITIAL SETTING OF THE IO CHANNELS
C  ALL INPUT CHANNELS ARE SET TO 5 ( STANDARD INPUT ) FOR VEGAS AND INPLOT
C  ALL OUTPUT CHANNELS ARE SET TO 6 ( STANDARD OUTPUT ) FOR VEGAS AND INPLOT
C  INPUT AND OUTPUT CHANNELS FOR SAVE AND RESTR ARE SET TO 7
C  INPUT AND OUTPUT CHANNELS FOR SAVE2 AND RESTR2 ARE SET TO 9
C  INPUT AND OUTPUT CHANNELS FOR VGSAVE AND VGRSTR ARE SET TO 8
C  INPUT AND OUTPUT CHANNELS FOR THE RANDOM NUMBER GENERATOR STATE IS SET TO 17
C  AN ADDITIONAL IO CHANNEL IS SET TO 18 ( SO THIS CHANNEL IS RESERVED FOR VEGAS
C                                         ALSO )
C
C  AUTHOR : S. DE JONG
C
      COMMON/VGASIO/NINP,NOUTP
      COMMON/VGPLIO/NIN,NOUT
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
      LOGICAL FIRST
C
      DATA FIRST /.TRUE./
C
      IF(FIRST)THEN
         FIRST=.FALSE.
         NINP =5
         NIN  =5
         NOUTP=6
         NOUT =6
         LUN1 =7
         LUN2 =9
         LUN3 =8
         LUN4 =17
         LUN5 =18
      ENDIF
C
      RETURN
C
      END
      SUBROUTINE VEGAS(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C
C  SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG*N
C    - BY G.P.  LEPAGE  SEPT 1976/ (REV) APR 1978
C
C    -FTN5 VERSION 21-8-1984
C    -HBOOK/HPLOT INTERFACE 6-1-1985
C
C  AUTHOR                                       : G. P. LEPAGE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FXN
      REAL ZZF1,ZZW
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/VGASIO/NINP,NOUTP
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      COMMON/VGRES/S1,S2,S3,S4
      COMMON/FER/AVGI
      REAL S1,S2,S3,S4
C
C
      DATA XL,XU/10*0.,10*1./
      DATA NDMX/50/,ALPH/1.5/,ONE/1./,MDS/1/
C
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      IF(NPRN.GT.0)THEN
         IPR=0
      ELSE
         IPR=1
      ENDIF
      NDO=1
      DO 1 J=1,NDIM
         XI(1,J)=ONE
1     CONTINUE
C
      ENTRY VEGAS1(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...INITIALIZES CUMMULATIVE VARIABLES,BUT NOT GRID
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS1 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
C
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
C
      ENTRY VEGAS2(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C...NO INITIALIZATION
      CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS2 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      ND=NDMX
      NG=1
      IF(MDS.NE.0)THEN
         NG=(NCALL/2.)**(1./NDIM)
         MDS=1
         IF((2*NG-NDMX).GE.0)THEN
            MDS=-1
            NPG=NG/NDMX+1
            ND=NG/NPG
            NG=NPG*ND
         ENDIF
      ENDIF
C
      K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
         DX(J)=XU(J)-XL(J)
         XJAC=XJAC*DX(J)
3     CONTINUE
C
C  REBIN PRESERVING BIN DENSITY
C
      IF(ND.NE.NDO)THEN
         RC=NDO/XND
         DO 7 J=1,NDIM
            K=0
            XN=0.
            DR=XN
            I=K
4           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
5           IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 5
            DO 6  I=1,NDM
               XI(I,J)=XIN(I)
6           CONTINUE
            XI(ND,J)=ONE
7        CONTINUE
         NDO=ND
      ENDIF
C
       IF(NPRN.NE.0.AND.NPRN.NE.10)WRITE(NOUTP,200)NDIM,CALLS,IT,ITMX
     + ,ACC,MDS,ND
       IF(NPRN.EQ.10)WRITE(NOUTP,290)NDIM,CALLS,ITMX,ACC,MDS,ND
C
      ENTRY VEGAS3(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C     - MAIN INTEGRATION LOOP
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)'0VEGAS3 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
9     CONTINUE
      IT=IT+1
      TI=0.
      TSI=TI
C
      DO 10 J=1,NDIM
         KG(J)=1
         DO 10 I=1,ND
            D(I,J)=TI
            DI(I,J)=TI
10    CONTINUE
C
11    FB=0.
      F2B=FB
      K=0
C
12    CONTINUE
      K=K+1
      DO 121 J=1,NDIM
         QRAN(J)=VGRAN(0.0)
121   CONTINUE
      WGT=XJAC
      DO 15 J=1,NDIM
         XN=(KG(J)-QRAN(J))*DXG+ONE
         IA(J)=XN
         IAJ=IA(J)
         IAJ1=IAJ-1
         IF(IAJ.LE.1)THEN
            XO=XI(IAJ,J)
            RC=(XN-IAJ)*XO
         ELSE
            XO=XI(IAJ,J)-XI(IAJ1,J)
            RC=XI(IAJ1,J)+(XN-IAJ)*XO
         ENDIF
         X(J)=XL(J)+RC*DX(J)
         WGT=WGT*XO*XND
15    CONTINUE
C
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
C
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
         IAJ=IA(J)
         DI(IAJ,J)=DI(IAJ,J)+F/CALLS
         IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
16    CONTINUE
      IF(K.LT.NPG) GO TO 12
C
      F2B=F2B*NPG
      F2B=DSQRT(F2B)
      F2B=DABS((F2B-FB)*(F2B+FB))
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.LT.0)THEN
         DO 17 J=1,NDIM
            IAJ=IA(J)
            D(IAJ,J)=D(IAJ,J)+F2B
17       CONTINUE
      ENDIF
      K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
C
C  FINAL RESULTS FOR THIS ITERATION
C
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      IF(TSI .EQ. 0.)THEN
         WGT = 0.
      ELSE
         WGT=TI2/TSI
      ENDIF
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      IF(SWGT .EQ. 0.)THEN
         AVGI=TI
      ELSE
         AVGI=SI/SWGT
      ENDIF
      IF(SI2 .EQ. 0.)THEN
         SD=TSI
      ELSE
         SD=SWGT*IT/SI2
      ENDIF
      SCALLS=SCALLS+CALLS
      CHI2A=0.
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      IF(SD .NE. 0.)THEN
         SD=ONE/SD
         SD=DSQRT(SD)
      ELSE
         SD=TSI
      ENDIF
      IF(NPRN.NE.0)THEN
         TSI=DSQRT(TSI)
         IF(NPRN.NE.10)WRITE(NOUTP,201)IPR,IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.EQ.10)WRITE(NOUTP,203)IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.LT.0)THEN
            DO 20 J=1,NDIM
               WRITE(NOUTP,202)J
               WRITE(NOUTP,204)(XI(I,J),DI(I,J),D(I,J),I=1,ND)
20         CONTINUE
         ENDIF
      ENDIF
C
C   REFINE GRID
C
21    IF(SD .NE. 0.)THEN
         REL = DABS(SD/AVGI)
      ELSE
         REL = 0.
      ENDIF
      IF(REL.LE.DABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
C
      DO 23 J=1,NDIM
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2.
         DT(J)=D(1,J)
         DO 22 I=2,NDM
            D(I,J)=XO+XN
            XO=XN
            XN=D(I+1,J)
            D(I,J)=(D(I,J)+XN)/3.
            DT(J)=DT(J)+D(I,J)
22       CONTINUE
         D(ND,J)=(XN+XO)/2.
         DT(J)=DT(J)+D(ND,J)
23    CONTINUE
C
      DO 28 J=1,NDIM
         RC=0.
         DO 24 I=1,ND
            R(I)=0.
            IF(D(I,J).GT.0.)THEN
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
            ENDIF
            RC=RC+R(I)
24       CONTINUE
         RC=RC/XND
         K=0
         XN=0.
         DR=XN
         I=K
25       K=K+1
         DR=DR+R(K)
         XO=XN
         XN=XI(K,J)
26       IF(RC.GT.DR) GO TO 25
         I=I+1
         DR=DR-RC
         IF(DR .EQ. 0.)THEN
            XIN(I)=XN
         ELSE
            XIN(I)=XN-(XN-XO)*DR/R(K)
         ENDIF
         IF(I.LT.NDM) GO TO 26
         DO 27 I=1,NDM
            XI(I,J)=XIN(I)
27       CONTINUE
         XI(ND,J)=ONE
28    CONTINUE
C
      IF(IT.LT.ITMX.AND.DABS(ACC).LT.REL)GO TO 9
C
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
C
199   FORMAT(A)
200   FORMAT('0INPUT PARAMETERS FOR VEGAS   NDIM=',I3
     +,'  NCALL=',F8.0/28X,'  IT=',I5,'  ITMX =',I5/28X
     +,'  ACC=',D9.3/28X,'  MDS=',I3,'   ND=',I4//)
290    FORMAT('0VEGAS  NDIM=',I3,'  NCALL=',F8.0,'  ITMX =',I5
     + ,'  ACC=',D9.3,'  MDS=',I3,'   ND=',I4)
201   FORMAT(/I1,'INTEGRATION BY VEGAS'/'0ITERATION NO',I3,
     +'.   INTEGRAL =',D14.8/20X,'STD DEV  =',D10.4/
     +' ACCUMULATED RESULTS.   INTEGRAL =',D14.8/
     +24X,'STD DEV  =',D10.4 / 24X,'CHI**2 PER ITN   =',D10.4)
202   FORMAT('0DATA FOR AXIS',I2 / 7X,'X',7X,'  DELT I  ',
     +2X,' CONVCE    ',11X,'X',7X,'  DELT I  ',2X,' CONVCE     '
     +,11X,'X',7X,'  DELT I  ',2X,' CONVCE     '/)
204   FORMAT(1X,3D12.4,5X,3D12.4,5X,3D12.4)
203   FORMAT(1X,I3,D20.8,D12.4,D20.8,D12.4,D12.4)
C
      END
      double precision FUNCTION VGRAN(S)
        double precision dranf
C
C   PICK YOUR CHOICE FOR A RANDOM NUMBER GENERATOR !
C
C  AUTHOR : J. VERMASEREN
C
C
      VGRAN  = DRANF(S)
      RETURN
      END
      SUBROUTINE RANDAT
C
C  INITIALISES THE NUMBER NCALL TO 0 TO FLAG THE FIRST CALL
C  OF THE RANDOM NUMBER GENERATOR
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
C
      LOGICAL FIRST
C
C
      DATA FIRST /.TRUE./
C
      IF(FIRST)THEN
         FIRST=.FALSE.
         NCALL=0
      ENDIF
C
      RETURN
C
      END
      FUNCTION RANF(DUMMY)
C
C   RANDOM NUMBER FUNCTION TAKEN FROM KNUTH
C   (SEMINUMERICAL ALGORITHMS).
C   METHOD IS X(N)=MOD(X(N-55)-X(N-24),1/FMODUL)
C   NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.
C
C   RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.
C   IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.
C   IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE.
C
C  AUTHOR                     : J. VERMASEREN
C  EXTENSION TO START THROUGH : S. DE JONG
C
      PARAMETER (FMODUL=1.E-09)
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
C
C
      CALL RANDAT
C
      IF( NCALL.EQ.0 ) THEN
          CALL IN55 ( IA,234612947 )
          MCALL = 55
          NCALL = 1
      ENDIF
      IF ( MCALL.EQ.0 ) THEN
          CALL IRN55(IA)
          MCALL=55
      ENDIF
      RANF=IA(MCALL)*FMODUL
      MCALL=MCALL-1
      RETURN
      END
      DOUBLE PRECISION FUNCTION DRANF(DUMMY)
C
C  AUTHOR : J. VERMASEREN
C
      EXTERNAL RANF
C
C
      DRANF=DBLE(RANF(DUMMY))
      RETURN
      END
      SUBROUTINE IN55(IA,IX)
C
C  AUTHOR : J. VERMASEREN
C
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
C
      IA(55)=IX
      J=IX
      K=1
      DO 10 I=1,54
         II=MOD(21*I,55)
         IA(II)=K
         K=J-K
         IF(K.LT.0)K=K+MODULO
         J=IA(II)
10    CONTINUE
      DO 20 I=1,10
         CALL IRN55(IA)
20    CONTINUE
      RETURN
      END
      SUBROUTINE IRN55(IA)
C
C  AUTHOR : J. VERMASEREN
C
      PARAMETER (MODULO=1000000000)
      INTEGER IA(55)
C
      DO 10 I=1,24
         J=IA(I)-IA(I+31)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
10    CONTINUE
      DO 20 I=25,55
         J=IA(I)-IA(I-24)
         IF(J.LT.0)J=J+MODULO
         IA(I)=J
20    CONTINUE
      RETURN
      END
      SUBROUTINE SIRN55
C
C  THIS ROUTINE SAVES THE STATE OF THE RANDOM NUMBER GENERATOR RANF
C  FROM LOGICAL UNIT LUN4
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL RANDAT
      CALL VGDAT
C
      WRITE(LUN4,1)NCALL,MCALL,(IA(I),I=1,55)
C
1     FORMAT(2I20,(5I16))
C
      RETURN
      END
      SUBROUTINE RIRN55
C
C  THIS ROUTINE READS THE STATE OF THE RANDOM NUMBER GENERATOR RANF
C  FROM LOGICAL UNIT LUN4
C
C  AUTHOR : S. DE JONG
C
      COMMON /CIRN55/NCALL,MCALL,IA(55)
      INTEGER IA
      COMMON/VGSAV/LUN1,LUN2,LUN3,LUN4,LUN5
C
C
      CALL RANDAT
      CALL VGDAT
C
      READ(LUN4,1)NCALL,MCALL,(IA(I),I=1,55)
C
1     FORMAT(2I20,(5I16))
C
      RETURN
      END










