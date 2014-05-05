!дома ночь 7-8 12
module generator    !модуль генератора
 integer(4) n1, n2, n3 !случаные числа
 integer(4) m1,a1,b1,a2,b2,m3,m2 !параметры генератора
 integer(4) max1, max2,max3 !максимумы генераторов
 real(8) randmass(500) !массив случайных чисел
 real(8) outrand
contains
subroutine randomn
!use generator
    implicit none

integer(4) i
real(8) sseed1,sseed2,sseed3
integer(4) seed(20)
character(20) a,b,c
real(8) n11,n21,n31
call date_and_time(a,b,c,seed)
sseed1=0.0
sseed2=0.0
sseed3=0.0
!print *,seed
do i=1,20,1 !здесь шарной бред
 sseed1=sseed1+seed(i)
 if (mod(i,2)==0) then
  sseed2=sseed2+seed(i)
 else
  sseed2=sseed2/2.0
 endif
 if (mod(i,3)==0) then
  sseed3=sseed3+seed(i)
 else
  sseed3=sseed3-seed(i)/2.2
 endif
enddo

!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
print *,'- Random seeds -', n11,n21,n31
!call random_seed(int(sseed))
!call random_number(n11)
!call random_number(n21)
!call random_number(n31)
m1=2147483563
a1=40014
b1=12345
m2=2147483399
a2=40692
b2=54321
max3=2147483647
max1=m1-1
max2=m2-1
n1=ceiling((n11*max1)-1)
n2=ceiling((n21*max2)-1)
n3=ceiling((n31*max3)-1)
!   open(15,file='rand.txt')
 do i=1,500
 call randstart()
 randmass(i)=outrand
 enddo
print *,'--- Generator Start ---'
end subroutine
!*******************************************************
!Рандомная функция начало
subroutine  randstart()
!use generator
    implicit none

n1=abs(mod(a1*n1+b1,m1))
n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
if (float(n3)/float(max3)<0.5) then
 outrand=float(n1)/float(max1)
else
 outrand=float(n2)/float(max2)
endif
end subroutine
!*******************************************************
!Рандомная функция
function getrand()
!use generator
    implicit none
real(8) getrand
integer(4) repick
!n1=abs(mod(a1*n1+b1,m1))
!n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
repick=ceiling((float(n3)/float(max3))*500)
getrand=randmass(repick)
call randstart()
randmass(repick)=outrand
end function

end module


program dcvmd
	use generator
	implicit none
	real(8) BoxL1,BoxL2
	real(8) BoxG1,BoxG2
	real(8) BoxH
	character(8) TempString,TempString1
	integer(4) TotalAtom
	integer(4) i,j,k
	integer(4) ii,jj,kk
	integer(4) SubNum
	character(5),allocatable:: SubName(:)
	integer(4),allocatable:: NLiq(:)
	integer(4),allocatable:: SubAtomNum(:)
	character(5),allocatable:: SubAtomName(:,:,:)
	real(8),allocatable:: SubAtomX(:,:,:),SubAtomY(:,:,:),SubAtomZ(:,:,:)
	real(8),allocatable:: SubAtomVX(:,:,:),SubAtomVY(:,:,:),SubAtomVZ(:,:,:)
	integer(4) MemAtom
	real(8),allocatable:: MemX(:),MemY(:),MemZ(:)
	character(5), allocatable:: MemName(:)
	integer(4) TempInt,TempInt1,TempInt2
	real(8) TempR1,TempR2,TempR3
	real(8) Sigma
	real(8) MemDelta
	integer(4) Mem1HW
	integer(4) Mem1Len
	real(8) RoL
	integer(4) TotalMol
	
	real(8),allocatable:: MolX(:,:),MolY(:,:),MolZ(:,:)
	integer(4),allocatable:: InOutV1(:,:),InOutV2(:,:)
	integer(4) MemType
	character(20),allocatable:: SubFile(:),SubFileINP(:)
	real(8),allocatable:: SubSigma(:,:)
	real(8),allocatable:: SubEps(:,:)
	real(8),allocatable:: SubMass(:,:)
	real(8),allocatable:: SubCha(:,:)
	
	real(8) SumX,SumY,SumZ,SumMass
	real(8),allocatable:: Vol1X(:,:,:),Vol1Y(:,:,:),Vol1Z(:,:,:)
	real(8),allocatable:: Vol2X(:,:,:),Vol2Y(:,:,:),Vol2Z(:,:,:)
	real(8),allocatable:: Vol3X(:,:,:),Vol3Y(:,:,:),Vol3Z(:,:,:)
	real(8),allocatable:: Vol1VX(:,:,:),Vol1VY(:,:,:),Vol1VZ(:,:,:)
	real(8),allocatable:: Vol2VX(:,:,:),Vol2VY(:,:,:),Vol2VZ(:,:,:)
	real(8),allocatable:: Vol3VX(:,:,:),Vol3VY(:,:,:),Vol3VZ(:,:,:)
	integer(4),allocatable:: SumInV1(:),SumInV2(:),SumOut(:)
	real(8),allocatable:: CenAtomX(:,:),CenAtomY(:,:),CenAtomZ(:,:)
	real(8),allocatable:: CheckX(:),CheckY(:),CheckZ(:)
	real(8) RandX,RandY,RandZ
	
	integer(4) NStep, step
	real(8) Vol1Ak,Vol2Ak
	integer(4) ch, CheckType
	real(8) DeltaEn,drx,dry,drz,rx,ry,rz,r,LJ
	real(8),allocatable:: MixSig(:,:,:,:),MixEps(:,:,:,:),MixCha(:,:,:,:)
	real(8) BoxVol,Temp,Prob
	real(8),allocatable:: AkL(:),AkG(:)
	real(8),allocatable:: Vol1MX(:,:),Vol1MY(:,:),Vol1MZ(:,:)
	real(8),allocatable:: Vol2MX(:,:),Vol2MY(:,:),Vol2MZ(:,:)
	real(8),allocatable:: Vol3MX(:,:),Vol3MY(:,:),Vol3MZ(:,:)
	
	integer(4) StepType,CheckMol
	real(8) rmx,rmy,rmz
	real(8) TotBoxLen
	integer(4) Ntot
	real(8),allocatable:: frac(:)
	character(20) time1,time2,time3
	integer(4) seed(20)
	integer(4) rseed
	integer(4),allocatable:: NMolLiq(:)
	real(8),allocatable:: ToRoLIq(:)
	character(100) LongTempString
	integer(4) TopInt
	real(8) CurTime
	integer(4) First
	real(8) DTime
	real(8) Mem1B,Mem1E,Mem2B,Mem2E
	
	call randomn
	
	open(7,file='test.temp')
		read(7,'(f20.10)') BoxL1
		read(7,'(f20.10)') BoxL2
		read(7,'(f20.10)') BoxG1
		read(7,'(f20.10)') BoxG2
		read(7,'(f20.10)') BoxH
		read(7,'(i6)') SubNum
		read(7,'(i6)') TotalMol
		allocate(SubName(SubNum+1))
		allocate(NLiq(SubNum+1))
		allocate(SubAtomNum(SubNum+1))
		allocate(SubFile(SubNum+1))
		allocate(SubFileINP(SubNum+1))
		allocate(frac(SubNum+1))
		allocate(AkG(SubNum+1))
		allocate(AkL(SubNum+1))
		allocate(ToRoLiq(SubNum+1))
		do i=1,SubNum
			read(7,'(a5,2i6)') SubName(i),NLiq(i),SubAtomNum(i)
		enddo
		!print *, SubName(i),NLiq(i),SubAtomNum(i)
		!pause
		read(7,'(i6)') MemAtom
		read(7,'(f20.10)') TotBoxLen
		read(7,'(f20.10)') Mem1B
		read(7,'(f20.10)') Mem1E
		read(7,'(f20.10)') Mem2B
		read(7,'(f20.10)') Mem2E
		read(7,'(f20.10)') CurTime
		read(7,'(i6)') First
		TotalMol=20000
		allocate(SubAtomName(SubNum+1,TotalMol,30))
		allocate(SubAtomX(SubNum+1,TotalMol,30))
		allocate(SubAtomY(SubNum+1,TotalMol,30))
		allocate(SubAtomZ(SubNum+1,TotalMol,30))
		allocate(SubAtomVX(SubNum+1,TotalMol,30))
		allocate(SubAtomVY(SubNum+1,TotalMol,30))
		allocate(SubAtomVZ(SubNum+1,TotalMol,30))
		allocate(MemName(MemAtom))
		allocate(MemX(MemAtom))
		allocate(MemY(MemAtom))
		allocate(MemZ(MemAtom))
	close(7)
	print *, 'Temp FILE ReAD DONE'
	
	open(9,file='main.in')
		read(9,*) TempString
		read(9,*) SubNum
		do i=1,SubNum
			read(9,*) SubFile(i),frac(i), AkL(i),AkG(i),ToRoLIq(i)
			write(SubFileINP(i),'(2a)') trim(adjustl(SubFile(i))),'.inp'
		enddo
		read(9,*) TempString
		read(9,*) MemType
		if (MemType==1) then
		read(9,*) TempString
			read(9,*) Sigma
			read(9,*) TempString
			read(9,*) MemDelta
			read(9,*) TempString
			read(9,*) Mem1HW
			read(9,*) TempString
			read(9,*) Mem1Len
			read(9,*) TempString
			read(9,*) RoL
		endif
		read(9,*) TempString
		read(9,*) NStep
		print *, TempString, NStep
		read(9,*) TempString
		read(9,*) Temp
		read(9,*) TempString
		read(9,*) TopInt
		read(9,*) TempString
		read(9,*) DTime
	close(9)
	print *, 'MEMIN File read DONE'
		
	allocate(SubSigma(SubNum+1,30))
	allocate(SubEps(SubNum+1,30))
	allocate(SubMass(SubNum+1,30))
	allocate(SubCha(SubNum+1,30))
	do i=1,SubNum
		print *,SubFile(i),SubAtomNum(i)
		open(41,file=SubFileINP(i))
		do j=1,SubAtomNum(i)
			read(41,*) SubSigma(i,j),SubEps(i,j),SubMass(i,j),SubCha(i,j)
			print *, i,j,SubSigma(i,j),SubEps(i,j),SubMass(i,j),SubCha(i,j)
		enddo
		close(41)
	enddo
	print *, ' Parameters file read DONE '
	!pause
		
	open(8,file='test.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, TempString,SubAtomName(i,j,k),&
					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
					print *,i,j,k, TempInt, TempString,SubAtomName(i,j,k),SubAtomName(i,1,k)
!					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
!					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
				enddo
			enddo
			!pause
		enddo
		do i=1,MemAtom
			read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt,TempString,MemName(i),&
			&TempInt1,MemX(i),MemY(i),MemZ(i),&
			&TempR1,TempR2,TempR3
		enddo
	close(8)
	print *, ' GRO FILE READ DONE '
	allocate(MolX(SubNum+1,TotalAtom))
	allocate(MolY(SubNum+1,TotalAtom))
	allocate(MolZ(SubNUm+1,TotalAtom))
	do i=1,SubNum
		do j=1,NLiq(i)
			SumX=0.0
			SumY=0.0
			SumZ=0.0
			SumMass=0.0
			do k=1,SubAtomNum(i)
				SumX=SumX+SubAtomX(i,j,k)*SubMass(i,k)
				SumY=SumY+SubAtomY(i,j,k)*SubMass(i,k)
				SumZ=SumZ+SubAtomZ(i,j,k)*SubMass(i,k)
				SumMass=SumMass+SubMass(i,k)
			enddo
			MolX(i,j)=SumX/SumMass
			MolY(i,j)=SumY/SumMass
			MolZ(i,j)=SumZ/SumMass
		enddo
	enddo
	
	print *, 'Molecule centers calculating DONE'
	allocate(SumInV1(SubNum+1))
	allocate(SumInV2(SubNum+1))
	allocate(SumOut(SubNum+1))
	allocate(InOutV1(SubNum+1,TotalAtom))
	allocate(InOutV2(SubNum+1,TotalAtom))
	do i=1,SubNum
		SumInV1(i)=0
		SumInV2(i)=0
		SumOut(i)=0
		do j=1,NLiq(i)
			InOutV1(i,j)=0
			InOutV2(i,j)=0
			if ((MolZ(i,j)>BoxL1).and.(MolZ(i,j)<BoxL2)) then
				InOutV1(i,j)=1
				SumInV1(i)=SumInV1(i)+1
			endif
			if ((MolZ(i,j)>BoxG1).and.(MolZ(i,j)<BoxG2)) then
				InOutV2(i,j)=1
				SumInV2=SumInV2+1
			endif
		enddo
		SumOut(i)=NLiq(i)-SumInV1(i)-SumInV2(i)
	enddo
	
	print *, 'Summing numbers of molecules DONE'
	allocate(Vol1X(SubNum+1,10000,10))
	allocate(Vol1Y(SubNum+1,10000,10))
	allocate(Vol1Z(subNum+1,10000,10))
	allocate(Vol1VX(SubNum+1,10000,10))
	allocate(Vol1VY(SubNum+1,10000,10))
	allocate(Vol1VZ(subNum+1,10000,10))
	
	allocate(Vol2X(SubNum+1,10000,10))
	allocate(Vol2Y(SubNum+1,10000,10))
	allocate(Vol2Z(SubNUm+1,10000,10))
	allocate(Vol2VX(SubNum+1,10000,10))
	allocate(Vol2VY(SubNum+1,10000,10))
	allocate(Vol2VZ(SubNum+1,10000,10))
	allocate(Vol3X(SubNum+1,10000,10))
	allocate(Vol3Y(SubNum+1,10000,10))
	allocate(Vol3Z(SubNum+1,10000,10))
	allocate(Vol3VX(SubNum+1,10000,10))
	allocate(Vol3VY(SubNum+1,10000,10))
	allocate(Vol3VZ(SubNum+1,10000,10))
	
	allocate(Vol1MX(SubNum+1,10000))
	allocate(Vol1MY(SubNum+1,10000))
	allocate(Vol1MZ(SubNum+1,10000))
	allocate(Vol2MX(SubNum+1,10000))
	allocate(Vol2MY(SubNum+1,10000))
	allocate(Vol2MZ(SubNum+1,10000))
	allocate(Vol3MX(SubNum+1,10000))
	allocate(Vol3MY(SubNum+1,10000))
	allocate(Vol3MZ(SubNum+1,10000))
	
	print*, 'Mega Allocating DONE'
	do i=1,SubNum
		TempInt=1
		TempInt1=1
		TempInt2=1
		do j=1,NLiq(i)
			if (InOutV1(i,j)==1) then
				Vol1MX(i,TempInt)=MolX(i,j)
				Vol1MY(i,TempInt)=MolY(i,j)
				Vol1MZ(i,TempInt)=MolZ(i,j)
				do k=1,SubAtomNum(i)
					Vol1X(i,TempInt,k)=SubAtomX(i,j,k)
					Vol1Y(i,TempInt,k)=SubAtomY(i,j,k)
					Vol1Z(i,TempInt,k)=SubAtomZ(i,j,k)
					Vol1VX(i,TempInt,k)=SubAtomVX(i,j,k)
					Vol1VY(i,TempInt,k)=SubAtomVY(i,j,k)
					Vol1VZ(i,TempInt,k)=SubAtomVZ(i,j,k)
				enddo
				TempInt=TempInt+1
			elseif (InOutV2(i,j)==1) then
				Vol2MX(i,TempInt1)=MolX(i,j)
				Vol2MY(i,TempInt1)=MolY(i,j)
				Vol2MZ(i,TempInt1)=MolZ(i,j)
				do k=1,SubAtomNum(i)
					Vol2X(i,TempInt1,k)=SubAtomX(i,j,k)
					Vol2Y(i,TempInt1,k)=SubAtomY(i,j,k)
					Vol2Z(i,TempInt1,k)=SubAtomZ(i,j,k)
					Vol2VX(i,TempInt1,k)=SubAtomVX(i,j,k)
					Vol2VY(i,TempInt1,k)=SubAtomVY(i,j,k)
					Vol2VZ(i,TempInt1,k)=SubAtomVZ(i,j,k)
				enddo
				TempInt1=TempInt1+1
			else
					Vol3MX(i,TempInt2)=MolX(i,j)
					Vol3MY(i,TempInt2)=MolY(i,j)
					Vol3MZ(i,TempInt2)=MolZ(i,j)
				do k=1,SubAtomNum(i)
					Vol3X(i,TempInt2,k)=SubAtomX(i,j,k)
					Vol3Y(i,TempInt2,k)=SubAtomY(i,j,k)
					Vol3Z(i,TempInt2,k)=SubAtomZ(i,j,k)
					Vol3VX(i,TempInt2,k)=SubAtomVX(i,j,k)
					Vol3VY(i,TempInt2,k)=SubAtomVY(i,j,k)
					Vol3VZ(i,TempInt2,k)=SubAtomVZ(i,j,k)
				enddo
				TempInt2=TempInt2+1
			endif
		enddo
	enddo
	
	print *, 'Rewriting coordinates by box parts'
	call date_and_time(time1,time2,time3,seed)
	rseed=ceiling(abs(cos(float(seed(8))/1000.0))*20000.0)
	print *,' Random integer ', rseed
	call srand(rseed)
	
	!centrating molecules
	allocate(CenAtomX(SubNum+1,30))
	allocate(CenAtomY(SubNum+1,30))
	allocate(CenAtomZ(SubNum+1,30))
	
	do i=1,SubNum
		do j=1,SubAtomNum(i)
			CenAtomX(i,j)=SubAtomX(i,1,j)-MolX(i,1)	!хотя бы одна молеула быть должна
			CenAtomY(i,j)=SubAtomY(i,1,j)-MolY(i,1)
			CenAtomZ(i,j)=SubAtomZ(i,1,j)-MolZ(i,1)
			print *,i,j, CenAtomX(i,j),CenAtomY(i,j),CenAtomZ(i,j)
		enddo
	enddo
	print *, ' To Center of Mass DONE'
	!pause
	allocate(CheckX(30))
	allocate(CheckY(30))
	allocate(CheckZ(30))
	
	allocate(MixSig(SubNum+1,SubNum+1,30,30))
	allocate(MixEps(SubNum+1,SubNum+1,30,30))
	allocate(MixCha(SubNum+1,SubNum+1,30,30))
	do i=1,SubNum
		do j=1,SubNum
			do k=1,SubAtomNum(i)
				do ch=1,SubAtomNum(j)
					MixSig(i,j,k,ch)=(SubSigma(i,k)+SubSigma(j,ch))/2.0
					MixEps(i,j,k,ch)=sqrt(SubEps(i,k)*SubEps(j,ch))
					MixCha(i,j,k,ch)=SubCha(i,k)*SubCha(j,ch)
					print *, i,k,k,ch, MixSig(i,j,k,ch),MixEps(i,j,k,ch),MixCha(i,j,k,ch)
				enddo
			enddo
		enddo
	enddo
	print *, ' CAlculating Mixrule DONE'
	!Vol1 create
	!pause
	
	BoxVol=BoxH*BoxH*BoxH
	!print*, log(AkL(1)*BoxVol/float(SumInV1(1)+1))
	
	do step=1,NStep
	StepType=ceiling(getrand()*4.0)
	StepType=1

	if (StepType==1) then
		RandX=getrand()*BoxH
		RandY=getrand()*BoxH
		RandZ=getrand()*(BoxH-0.6)+BoxL1+0.3
		CheckType=ceiling(getrand()*float(SubNum))
		if (SumInV1(CheckType)/BoxVol*Sigma*Sigma*Sigma<ToRoLIq(CheckType)) then
			do i=1,SubAtomNum(CheckType)
				CheckX(i)=RandX+CenAtomX(CheckType,i)
				CheckY(i)=RandY+CenAtomY(CheckType,i)
				CheckZ(i)=RandZ+CenAtomZ(CheckType,i)
			enddo
			DeltaEn=0.0
			do i=1,SubNum
				do j=1,SumInV1(i)
!					rmx=abs(RandX-Vol1MX(i,j))
!					rmy=abs(RandY-Vol1MY(i,j))
!					rmz=abs(RandZ-Vol1MZ(i,j))
!					if (rmx>BoxH/2.0) then
!						drx=BoxH/2.0
!					else
!						drx=0.0
!					endif
!					if (rmy>BoxH/2.0) then
!						dry=BoxH/2.0
!					else 
!						dry=0.0
!					endif
!					if (rmz>BoxH/2.0) then
!						drz=BoxH/2.0
!					else 
!						drz=0.0
!					endif
!					rmx=rmx-drx
!					rmy=rmy-dry
!					rmz=rmz-drz
!					if (sqrt(rmx*rmx+rmy*rmy+rmz*rmz)<BoxH/2.0) then
						do k=1,SubAtomNum(i)
							do ch=1,SubAtomNum(CheckType)
								rx=abs(CheckX(ch)-Vol1X(i,j,k)) !-drx
								if (rx>BoxH/2.0) then
									rx=BoxH-rx
								endif
								ry=abs(CheckY(ch)-Vol1Y(i,j,k)) !-dry
								if (ry>BoxH/2.0) then
									ry=BoxH-ry
								endif
								rz=abs(CheckZ(ch)-Vol1Z(i,j,k)) !-drz
								!if (rz>BoxH/2.0) then
								!	rz=BoxH-rz
								!endif
								r=sqrt(rx*rx+ry*ry+rz*rz)
								rz=MixSig(i,CheckType,k,ch)/r
								rz=rz*rz
								rz=rz*rz*rz
								if(r<0.23) then
									LJ=99999999999.0
								else
								!	LJ=4.0*MixEps(i,CheckType,k,ch)*(rz*rz-rz)	!&
								!&+MixCha(i,CheckType,k,ch)/r*167100.96
								endif
								!print *,r, LJ
								!pause
								DeltaEn=DeltaEn+LJ
							enddo
						enddo
!					endif
				enddo
			enddo
		
			Prob=exp(-DeltaEn/Temp)   !!*AkL(CheckType)     !log(AkL(CheckType)*BoxVol/float(SumInV1(CheckType)+1)))
			!Prob=2.0
			if (DeltaEn<AkL(1)) then
				print *,step, ' Molecule type ', CheckType, ' number ', SumInV1(CheckType)+1, ' add to VOL 1'
				!print*,'Prob ', AkL(1),'Delta En ', DeltaEn,log(AkL(CheckType)*BoxVol/(SumInV1(CheckType)+1))
				SumInV1(CheckType)=SumInV1(CheckType)+1
				Vol1MX(CheckType,SumInV1(CheckType))=RandX
				Vol1MY(CheckType,SumInV1(CheckType))=RandY
				Vol1MZ(CheckType,SumInV1(CheckType))=RandZ
				do i=1,SubAtomNum(CheckType)
					Vol1X(CheckType,SumInV1(CheckType),i)=CheckX(i)
					Vol1Y(CheckType,SumInV1(CheckType),i)=CheckY(i)
					Vol1Z(CheckType,SumInV1(CheckType),i)=CheckZ(i)
					Vol1VX(CheckType,SumInV1(CheckType),i)=0.0
					Vol1VY(CheckType,SumInV1(CheckType),i)=0.0
					Vol1VZ(CheckType,SumInV1(CheckType),i)=0.0
				enddo
			endif
		endif
	endif
	
	if (StepType==2) then
		RandX=getrand()*BoxH
		RandY=getrand()*BoxH
		RandZ=getrand()*BoxH+BoxG1
		CheckType=ceiling(rand()*float(SubNum))
		do i=1,SubAtomNum(CheckType)
			CheckX(i)=RandX+CenAtomX(CheckType,i)
			CheckY(i)=RandY+CenAtomY(CheckType,i)
			CheckZ(i)=RandZ+CenAtomZ(CheckType,i)
		enddo
		DeltaEn=0.0
		do i=1,SubNum
			do j=1,SumInV2(i)
				rmx=abs(RandX-Vol2MX(i,j))
				rmy=abs(RandY-Vol2MY(i,j))
				rmz=abs(RandZ-Vol2MZ(i,j))
				if (rmx>BoxH/2.0) then
					drx=BoxH/2.0
				else
					drx=0.0
				endif
				if (rmy>BoxH/2.0) then
					dry=BoxH/2.0
				else 
					dry=0.0
				endif
				if (rmz>BoxH*1.5) then
					drz=BoxH*1.5
				else 
					drz=0.0
				endif
				rmx=rmx-drx
				rmy=rmy-dry
				rmz=rmz-drz
				if (sqrt(rmx*rmx+rmy*rmy+rmz*rmz)<BoxH/2.0) then
					do k=1,SubAtomNum(i)
						do ch=1,SubAtomNum(CheckType)
							rx=abs(CheckX(ch)-Vol2X(i,j,k))-drx
							ry=abs(CheckY(ch)-Vol2Y(i,j,k))-dry
							rz=abs(CheckZ(ch)-Vol2Z(i,j,k))-drz
							r=sqrt(rx*rx+ry*ry+rz*rz)
							rz=MixSig(i,CheckType,k,ch)/r
							rz=rz*rz
							rz=rz*rz*rz
							LJ=4.0*MixEps(i,CheckType,k,ch)*(rz*rz-rz)&
							&+MixCha(i,CheckType,k,ch)/r*167100.96
							DeltaEn=DeltaEn+LJ
						enddo
					enddo
				endif
			enddo
		enddo
		Prob=exp(-DeltaEn/Temp+log(AkG(CheckType)*BoxVol*3.0/(SumInV2(CheckType)+1)))
		if (Prob>rand()) then
			print *, 'Molecule type ', CheckType, ' number ', SumInV2(CheckType)+1, ' add to VOL2 '
			SumInV2(CheckType)=SumInV2(CheckType)+1
			Vol2MX(CheckType,SumInV1(CheckType)+1)=RandX
			Vol2MY(CheckType,SumInV1(CheckType)+1)=RandY
			Vol2MZ(CheckType,SumInV1(CheckType)+1)=RandZ
			do i=1,SubAtomNum(CheckType)
				Vol2X(CheckType,SumInV1(CheckType)+1,i)=CheckX(i)
				Vol2Y(CheckType,SumInV1(CheckType)+1,i)=CheckY(i)
				Vol2Z(CheckType,SumInV1(CheckType)+1,i)=CheckZ(i)
				Vol2VX(CheckType,SumInV1(CheckType)+1,i)=0.0
				Vol2VY(CheckType,SumInV1(CheckType)+1,i)=0.0
				Vol2VZ(CheckType,SumInV1(CheckType)+1,i)=0.0
			enddo
		endif
	endif
	
	if (StepType==3) then 
		CheckType=ceiling(rand()*float(SubNum))
		CheckMol=ceiling(rand()*float(SumInV1(CheckType)))
		DeltaEn=0.0
		do i=1,SubNum
			do j=1,SubAtomNum(i)
				if ((CheckType/=i).and.(CheckMol/=j)) then
					rmx=abs(Vol1MX(i,j)-Vol1MX(CheckType,CheckMol))
					rmy=abs(Vol1MY(i,j)-Vol1MY(CheckType,CheckMol))
					rmz=abs(Vol1MZ(i,j)-Vol1MZ(CheckType,CheckMol))
					if (rmx>BoxH/2.0) then
						drx=BoxH/2.0
					else
						drx=0.0
					endif
					if (rmy>BoxH/2.0) then
						dry=BoxH/2.0
					else 
						dry=0.0
					endif
					if (rmz>BoxH/2.0) then
						dry=BoxH/2.0
					else 
						dry=0.0
					endif
					rmx=rmx-drx
					rmy=rmy-dry
					rmz=rmz-drz
					if (sqrt(rmx*rmx+rmy*rmy+rmz*rmz)>BoxH/2.0) then
						do k=1,SubAtomNum(i)
							do ch=1,SubAtomNum(CheckType)
								rx=abs(Vol1X(i,j,k)-Vol1X(CheckType,CheckMol,ch))-drx
								ry=abs(Vol1Y(i,j,k)-Vol1Y(CheckType,CheckMol,ch))-dry
								rz=abs(Vol1Z(i,j,k)-Vol1Z(CheckType,CheckMol,ch))-drz
								r=sqrt(rx*rx+ry*ry+rz*rz)
								rz=MixSig(CheckType,i,j,ch)/r
								rz=rz*rz
								rz=rz*rz*rz
								LJ=4.0*MixEps(CheckType,i,j,ch)*(rz*rz-rz)&
								&+MixCha(CheckType,i,j,ch)/r*167100.96
								DeltaEn=DeltaEn+LJ
							enddo
						enddo
					endif
				endif
			enddo
		enddo
		Prob=exp(-DeltaEn/Temp+log(SumInV1(CheckType)/AkL(CheckType)/BoxVol))
		if (Prob<rand()) then
			print *, ' Molecule tupe ', CheckType, ' number ', CheckMol, ' delete from VOL1'
			do i=CheckMol,SumInV1(CheckType)
				do j=1,SubAtomNum(CheckType)
					Vol1X(CheckType,i,j)=Vol1X(CheckType,i+1,j)
					Vol1Y(CheckType,i,j)=Vol1Y(CheckType,i+1,j)
					Vol1Z(CheckType,i,j)=Vol1Z(CheckType,i+1,j)
					Vol1VX(CheckType,i,j)=Vol1VX(CheckType,i+1,j)
					Vol1VY(CheckType,i,j)=Vol1VY(CheckType,i+1,j)
					Vol1VZ(CheckType,i,j)=Vol1VZ(CheckType,i+1,j)
				enddo
				Vol1MX(CheckType,i)=Vol1MX(CheckType,i+1)
				Vol1MY(CheckType,i)=Vol1MY(CheckType,i+1)
				Vol1MZ(CheckType,i)=Vol1MZ(CheckType,i+1)
			enddo
			SumInV1(CheckType)=SumInV1(CheckType)-1
		endif
	endif
	
	if (StepType==4) then 
		CheckType=ceiling(rand()*float(SubNum))
		CheckMol=ceiling(rand()*float(SumInV2(CheckType)))
		DeltaEn=0.0
		do i=1,SubNum
			do j=1,SubAtomNum(i)
				if ((CheckType/=i).and.(CheckMol/=j)) then
					rmx=abs(Vol2MX(i,j)-Vol2MX(CheckType,CheckMol))
					rmy=abs(Vol2MY(i,j)-Vol2MY(CheckType,CheckMol))
					rmz=abs(Vol2MZ(i,j)-Vol2MZ(CheckType,CheckMol))
					if (rmx>BoxH/2.0) then
						drx=BoxH/2.0
					else
						drx=0.0
					endif
					if (rmy>BoxH/2.0) then
						dry=BoxH/2.0
					else 
						dry=0.0
					endif
					if (rmz>BoxH*1.5) then
						dry=BoxH*1.5
					else 
						dry=0.0
					endif
					rmx=rmx-drx
					rmy=rmy-dry
					rmz=rmz-drz
					if (sqrt(rmx*rmx+rmy*rmy+rmz*rmz)>BoxH/2.0) then
						do k=1,SubAtomNum(i)
							do ch=1,SubAtomNum(CheckType)
								rx=abs(Vol2X(i,j,k)-Vol2X(CheckType,CheckMol,ch))-drx
								ry=abs(Vol2Y(i,j,k)-Vol2Y(CheckType,CheckMol,ch))-dry
								rz=abs(Vol2Z(i,j,k)-Vol2Z(CheckType,CheckMol,ch))-drz
								r=sqrt(rx*rx+ry*ry+rz*rz)
								rz=MixSig(CheckType,i,j,ch)/r
								rz=rz*rz
								rz=rz*rz*rz
								LJ=4.0*MixEps(CheckType,i,j,ch)*(rz*rz-rz)&
								&+MixCha(CheckType,i,j,ch)/r*167100.96
								DeltaEn=DeltaEn+LJ
							enddo
						enddo
					endif
				endif
			enddo
		enddo
		Prob=exp(-DeltaEn/Temp+log(SumInV2(CheckType)/AkG(CheckType)/(BoxVol*3.0)))
		if (Prob<rand()) then
		print *, ' Molecule tupe ', CheckType, ' number ', CheckMol, ' delete from VOL2'
			do i=CheckMol,SumInV2(CheckType)
				do j=1,SubAtomNum(CheckType)
					Vol2X(CheckType,i,j)=Vol2X(CheckType,i+1,j)
					Vol2Y(CheckType,i,j)=Vol2Y(CheckType,i+1,j)
					Vol2Z(CheckType,i,j)=Vol2Z(CheckType,i+1,j)
					Vol2VX(CheckType,i,j)=Vol2VX(CheckType,i+1,j)
					Vol2VY(CheckType,i,j)=Vol2VY(CheckType,i+1,j)
					Vol2VZ(CheckType,i,j)=Vol2VZ(CheckType,i+1,j)
				enddo
				Vol2MX(CheckType,i)=Vol2MX(CheckType,i+1)
				Vol2MY(CheckType,i)=Vol2MY(CheckType,i+1)
				Vol2MZ(CheckType,i)=Vol2MZ(CheckType,i+1)
			enddo
			SumInV2(CheckType)=SumInV2(CheckType)-1
		endif
	endif
	
	enddo
	
	Ntot=0
	do i=1,SubNum
		do j=1,SumOut(i)
			Ntot=Ntot+SubAtomNum(i)
		enddo
		do j=1,SumInV1(i)
			Ntot=Ntot+SubAtomNum(i)
		enddo
!		do j=1,SumInV2(i)
!			Ntot=Ntot+1
!		enddo
	enddo
	
	!groout
	open(31,file='testout.gro')
	write(31,'(a)') ' generated by dcv '
	write(31,'(i6)') Ntot+MemAtom
	TempInt=1
	do i=1,SubNum
		do j=1,SumOut(i)
			do k=1,SubAtomNum(i)
				write(31,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, SubName(i) , SubAtomName(i,1,k),k,&
				&Vol3X(i,j,k),Vol3Y(i,j,k),Vol3Z(i,j,k),Vol3VX(i,j,k),Vol3VY(i,j,k),Vol3VZ(i,j,k)
				TempInt=TempInt+1
			enddo
		enddo
		do j=1,SumInV1(i)
			do k=1,SubAtomNum(i)
				print * ,i,j,k,  TempInt, SubName(i) , SubAtomName(i,j,k), SubAtomName(i,1,k)
				write(31,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, SubName(i) , SubAtomName(i,1,k),k,&
				&Vol1X(i,j,k),Vol1Y(i,j,k),Vol1Z(i,j,k),Vol1VX(i,j,k),Vol1VY(i,j,k),Vol1VZ(i,j,k)
				TempInt=TempInt+1
			enddo
		enddo
	enddo
!	do i=1,SubNum
!		do j=1,SumInV1(i)
!			do k=1,SubAtomNum(i)
!				write(31,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, SubName(i) , SubAtomName(i,1,k),k,&
!				&Vol1X(i,j,k),Vol1Y(i,j,k),Vol1Z(i,j,k),Vol1VX(i,j,k),Vol1VY(i,j,k),Vol1VZ(i,j,k)
!				TempInt=TempInt+1
!			enddo
!		enddo
!	enddo
!	do i=1,SubNum
!		do j=1,SumInV2(i)
!			do k=1,SubAtomNum(i)
!				write(31,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, SubName(i) , SubAtomName(i,1,k),k,&
!				&Vol2X(i,j,k),Vol2Y(i,j,k),Vol2Z(i,j,k),Vol2VX(i,j,k),Vol2VY(i,j,k),Vol2VZ(i,j,k)
!				TempInt=TempInt+1
!			enddo
!		enddo
!	enddo
	do i=1,MemAtom
		write(31,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt,'MEM',MemName(i),&
		&1,MemX(i),MemY(i),MemZ(i),&
		&0.0,0.0,0.0
		TempInt=TempInt+1
	enddo
	write(31,'(3f20.5)') BoxH, BoxH, TotBoxLen
	close(31)
	
	TotalMol=0		!здесь все меняется ОЛОЛО
	allocate(NMolLiq(SubNum+1))
	TotalMol=0
	do i=1,SubNum
		NLiq(i)=0
		NMolLiq(i)=0
		do j=1,SumInV1(i)
			NMolLiq(i)=NMolLiq(i)+1
			do k=1,SubAtomNum(i)
				
				TotalMol=TotalMol+1
			enddo
			NLiq(i)=NLiq(i)+1
		enddo
		do j=1,SumOut(i)
			NMolLiq(i)=NMolLiq(i)+1
			do k=1,SubAtomNum(i)
				
				TotalMol=TotalMol+1
			enddo
			NLiq(i)=NLiq(i)+1
		enddo
	enddo
	
	!correcting test.temp
	open(33,file='test.temp')
		write(33,'(f20.10)') BoxL1
		write(33,'(f20.10)') BoxL2
		write(33,'(f20.10)') BoxG1
		write(33,'(f20.10)') BoxG2
		write(33,'(f20.10)') BoxH
		write(33,'(i6)') SubNum
		write(33,'(i6)') TotalMol
		do i=1,SubNum
			write(33,'(a5,i6,i6)') SubName(i),NLiq(i),SubAtomNum(i)
		enddo
		write(33,'(i6)') MemAtom
		write(33,'(f20.10)') TotBoxLen
		write(33,'(f20.10)') Mem1B
		write(33,'(f20.10)') Mem1E
		write(33,'(f20.10)') Mem2B
		write(33,'(f20.10)') Mem2E
		write(33,'(f20.10)') CurTime+DTime
		write(33,'(i6)') First
	close(33)
	!change acrivity
	do i=1,SubNum
		if (SumInV1(i)/BoxVol*Sigma*Sigma*Sigma<0.98*ToRoLIq(i)) then
			AkL(i)=AkL(i)*1.05
		elseif (SumInV1(i)/BoxVol*Sigma*Sigma*Sigma>1.1*ToRoLIq(i)) then
			AkL(i)=AkL(i)*0.98
		endif
!		if (AkL(i)>10000.0) then
!			AkL(i)=10000.0
!		endif
	enddo
	
	TotalMol=0
	do i=1,SubNum
		NLiq(i)=0
		NMolLiq(i)=0
		do j=1,SumInV1(i)
			NMolLiq(i)=NMolLiq(i)+1
			do k=1,SubAtomNum(i)
				NLiq(i)=NLiq(i)+1
				TotalMol=TotalMol+1
			enddo
			
		enddo
		do j=1,SumOut(i)
			NMolLiq(i)=NMolLiq(i)+1
			do k=1,SubAtomNum(i)
				NLiq(i)=NLiq(i)+1
				TotalMol=TotalMol+1
			enddo
			
		enddo
	enddo
	
	open(15,file='dens.out')
		do i=1,SubNum
			write(15,'(a,f20.10)') SubName(i), SumInV1(i)/BoxVol*Sigma*Sigma*Sigma
		enddo
		write(15,'(f20.10)') BoxH,BoxVol
	close(15)
	
	open(9,file='main.in')
		write(9,'(a)') 'Numbers of Substanses'
		write(9,'(i6)') SubNum
		do i=1,SubNum
			write(9,'(a,4f20.10)') trim(adjustl(SubFile(i))),frac(i), AkL(i),AkG(i),ToRoLIq(i)
		enddo
		write(9,'(a)') 'Initial Mem Type'
		write(9,'(i6)') MemType
		if (MemType==1) then
			write(9,'(a)') 'Sigma in GROMACS'
			write(9,'(f10.5)') Sigma
			write(9,'(a)') 'MEMcell len'
			write(9,'(f10.5)') MemDelta
			write(9,'(a)') 'MEM H/W'
			write(9,'(i6)') Mem1HW
			write(9,'(a)') 'MEM Len'
			write(9,'(i6)') Mem1Len
			write(9,'(a)') 'Liq density'
			write(9,'(f10.5)') RoL
		endif
		write(9,'(a)') 'NSteps'
		write(9,'(i10)') NStep
		write(9,'(a)') 'Temperature'
		write(9,'(f10.5)') Temp
		write(9,'(a)') 'Top String'
		write(9,'(i6)') TopInt
		write(9,'(a)') 'Delta Time'
		write(9,'(f20.10)') DTime
	close(9)
	
	!change top
	open(41,file='test.top')
	open(42,file='testnew.top')
	do i=1,TopInt-1
		read(41,'(a)') LongTempString
		write(42,'(a)') adjustr(trim(adjustl(LongTempString)))
	enddo
	do i=1,SubNum
		write(42,'(a10,i6)') SubName(i),NMolLiq(i)
	enddo
	if (MemType==1) then
		write(42,'(a10,i6)') 'MEM', Mem1HW*Mem1HW*Mem1Len*2    !MemAtom
	else
		write(42,'(a10,i6)') 'MEM', 2 
	endif
	close(42)
	close(41)
	
	!change index file
	open(45,file='indexnew.ndx')
	TempInt=1
	TempInt1=1
	write(45,'(a)') '[ System ]'
	do i=1,SubNum
		do j=1,NLiq(i)
			write(45,'(i8,a1,$)') TempInt,' '
			if (mod(TempInt1,10)==0) then
				write(45,'(a)') ' '
			endif
			TempInt=TempInt+1
			TempInt1=TempInt1+1
		enddo
	enddo
	do i=1,MemAtom
		write(45,'(i8,a1,$)') TempInt,' '
		if (mod(TempInt,10)==0) then
			write(45,'(a)') ' '
		endif
		TempInt=TempInt+1
		TempInt1=TempInt1+1
	enddo
	
	TempInt=1
	write(45,'(a)') ' '
	do i=1,SubNum
		TempInt1=1
		write(45,'(3a)') '[',trim(adjustl(SubName(i))),' ]'
		do j=1,NLiq(i)
			write(45,'(i8,a1,$)') TempInt,' '
			if (mod(TempInt1,10)==0) then
				write(45,'(a)') ' '
			endif
			TempInt=TempInt+1
			TempInt1=TempInt1+1
		enddo
		write(45,'(a)') ' '
	enddo
	write(45,'(a)') ' '
	write(45,'(a)') '[ MEM ]'
	TempInt1=1
	do i=1,MemAtom
		write(45,'(i8,a1,$)') TempInt,' '
		if (mod(TempInt,10)==0) then
			write(45,'(a)') ' '
		endif
		TempInt=TempInt+1
		TempInt1=TempInt1+1
	enddo
	write(45,'(a)') '  '
	
	close(45)
	
end program

