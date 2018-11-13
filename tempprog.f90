program tempprog
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
	integer(4),allocatable:: SubAtomIndex(:,:,:)
	integer(4),allocatable:: SubAtomIndex2(:,:,:)

	integer(4) MemAtom
	real(8),allocatable:: MemX(:),MemY(:),MemZ(:)
	character(5), allocatable:: MemName(:)
	integer(4), allocatable:: MemIndex(:)
	integer(4), allocatable:: MemIndex2(:)

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
	
	real(8),allocatable:: Vol1MolX(:,:),Vol1MolY(:,:),Vol1MolZ(:,:)
	real(8),allocatable:: Vol2MolX(:,:),Vol2MolY(:,:),Vol2MolZ(:,:)
	real(8) Mem1B,Mem1E,Mem2B,Mem2E
	real(8) CurTime
	integer(4) First
	real(8) DTime
	
	integer(4) Nslice
	real(8),allocatable:: MolNumSl(:,:)
	integer(4) hist
	
	real(8),allocatable:: InMem1(:),MemIn1(:),MemOut1(:)
	real(8),allocatable:: InMem2(:),MemIn2(:),MemOut2(:)
	real(8),allocatable:: Vol1MolNum(:),Vol2MolNum(:)
	character(5),allocatable:: SubMolName(:,:,:)
	
	integer(4) SampleNum

	integer(4) WallAtom
	integer(4),allocatable:: WallIndex(:)
	integer(4),allocatable:: WallIndex2(:)
	real(8),allocatable:: WallX(:)
	real(8),allocatable:: WallY(:)
	real(8),allocatable:: WallZ(:)
	
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
		read(7,'(i6)') MemAtom
		read(7,'(f20.10)') TotBoxLen
		read(7,'(f20.10)') Mem1B
		read(7,'(f20.10)') Mem1E
		read(7,'(f20.10)') Mem2B
		read(7,'(f20.10)') Mem2E
		read(7,'(f20.10)') CurTime
		print *,CurTime
		read(7,'(i6)') First
		allocate(SubAtomName(SubNum+1,TotalMol,30))
		allocate(SubMolName(SubNum+1,TotalMol,30))
		allocate(SubAtomX(SubNum+1,TotalMol,30))
		allocate(SubAtomY(SubNum+1,TotalMol,30))
		allocate(SubAtomZ(SubNum+1,TotalMol,30))
		allocate(SubAtomVX(SubNum+1,TotalMol,30))
		allocate(SubAtomVY(SubNum+1,TotalMol,30))
		allocate(SubAtomVZ(SubNum+1,TotalMol,30))

		allocate(SubAtomIndex(SubNum+1,TotalMol,30))
		allocate(SubAtomIndex2(SubNum+1,TotalMol,30))
		
		allocate(MemName(MemAtom))
		allocate(MemIndex(MemAtom))
		allocate(MemIndex2(MemAtom))
		allocate(MemX(MemAtom))
		allocate(MemY(MemAtom))
		allocate(MemZ(MemAtom))
	close(7)
	print *, 'MAIN IN FILE ReAD DONE'
	
	open(9,file='main.in')
		read(9,*) TempString
		read(9,*) SubNum
		do i=1,SubNum
			read(9,*) SubFile(i),frac(i), AkL(i),AkG(i),ToRoLIq(i)
			write(SubFileINP(i),'(2a)') trim(adjustl(SubFile(i))),'.inp'
		enddo
		read(9,*) TempString
		read(9,*) MemType
!		if (MemType==1) then
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
!		endif
		read(9,*) TempString
		read(9,*) NStep
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
	
	open(8,file='test.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		print *, TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
				print * ,i,j,k
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') SubAtomIndex2(i,j,k), TempString,SubAtomName(i,j,k),&
					&SubAtomIndex(i,j,k), SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
					&TempR1,TempR2,TempR3
!					print *,i,j,k, TempInt, TempString,SubAtomName(i,j,k),&
!					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
!					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
				enddo
			enddo
		enddo
		do i=1,MemAtom
			read(8,'(i5,2a5,i5,3f8.3,3f8.4)') MemIndex2(i),TempString,MemName(i),&
			&MemIndex(i),MemX(i),MemY(i),MemZ(i),&
			&TempR1,TempR2,TempR3
		enddo
		! reaad wall atoms
		if(MemType .eq. 5) then
			WallAtom=Mem1HW*Mem1HW*4
			allocate(WallIndex(WallAtom))
			allocate(WallIndex2(WallAtom))
			allocate(WallX(WallAtom))
			allocate(WallY(WallAtom))
			allocate(WallZ(WallAtom))
			do i=1,WallAtom
				read(8,'(i5,2a5,i5,3f8.3,3f8.4)') WallIndex2(i),TempString,TempString1,&
				&WallIndex(i),WallX(i),WallY(i),WallZ(i),&
				&TempR1,TempR2,TempR3
				!print *, TempInt1, TempInt1, WallX(i), WallY(i), WallZ(i)
			enddo
		endif

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
	
	!замена координат
	
	allocate(Vol1MolX(SubNum+1,TotalAtom))
	allocate(Vol1MolY(SubNum+1,TotalAtom))
	allocate(Vol1MolZ(SubNum+1,TotalAtom))
	allocate(Vol2MolX(SubNum+1,TotalAtom))
	allocate(Vol2MolY(SubNum+1,TotalAtom))
	allocate(Vol2MolZ(SubNum+1,TotalAtom))
	
	do i=1,SubNum
		do j=1,NLiq(i)
			Vol1MolX(i,j)=MolX(i,j)
			Vol1MolY(i,j)=MolY(i,j)
			Vol1MolZ(i,j)=MolZ(i,j)
		enddo
	enddo
	
	open(8,file='testout.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, SubMolName(i,j,k),SubAtomName(i,j,k),&
					&TempInt1,TempR1,TempR2,TempR3,&
					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
!					print *,i,j,k, TempInt, TempString,SubAtomName(i,j,k),&
!					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
!					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
				enddo
			enddo
		enddo
		do i=1,MemAtom
			read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt,TempString,MemName(i),&
			&TempInt1,MemX(i),MemY(i),MemZ(i),&
			&TempR1,TempR2,TempR3
		enddo
	close(8)
	print *, ' GRO FILE READ DONE '
	
	!
	TempInt=1
	open(38,file='testnew.gro')
		write(38,'(a)') TempString
		write(38,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					write(38,'(i5,2a5,i5,3f8.3,3f8.4)') SubAtomIndex2(i,j,k), SubMolName(i,j,k),SubAtomName(i,j,k),&
					&SubAtomIndex(i,j,k),SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
					&SubAtomVX(i,j,k),SubAtomVY(i,j,k),SubAtomVZ(i,j,k)
					TempInt=TempInt+1
				enddo
			enddo
		enddo
		do i=1,MemAtom
			write(38,'(i5,2a5,i5,3f8.3,3f8.4)') MemIndex2(i),TempString,MemName(i),&
			&MemIndex(i),MemX(i),MemY(i),MemZ(i),&
			&TempR1,TempR2,TempR3
			TempInt=TempInt+1
		enddo
		if(MemType .eq. 5) then
			do i=1,WallAtom
				write(38,'(i5,2a5,i5,3f8.3,3f8.4)') WallIndex(i),'WAL','W',&
				&WallIndex(i),WallX(i),WallY(i),WallZ(i),&
				&TempR1,TempR2,TempR3
				TempInt=TempInt+1
			enddo
		endif
		
		write(38,'(3f20.5)') BoxH, BoxH, TotBoxLen
	close(38)
end program

