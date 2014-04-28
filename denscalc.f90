program dcvcalc
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
	real(8),allocatable:: NearInVol(:),NearOutVol(:)
	
	integer(4) SampleNum
	
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
	
	open(8,file='test1.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, TempString,SubAtomName(i,j,k),&
					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
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
	
	allocate(NearInVol(SubNum+1))
	allocate(NearOutVol(SubNum+1))
	
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
	
	open(8,file='test2.gro')
		read(8,'(a)') TempString
		read(8,'(i6)') TotalAtom
		do i=1,SubNum
			do j=1,NLiq(i)
				do k=1,SubAtomNum(i)
					read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, TempString,SubAtomName(i,j,k),&
					&TempInt1,SubAtomX(i,j,k),SubAtomY(i,j,k),SubAtomZ(i,j,k),&
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
	
	do i=1,SubNum
		do j=1,NLiq(i)
			Vol2MolX(i,j)=MolX(i,j)
			Vol2MolY(i,j)=MolY(i,j)
			Vol2MolZ(i,j)=MolZ(i,j)
		enddo
	enddo
	
	
	
	Nslice=ceiling(TotBoxLen/0.1)
	open(56,file='nslice.temp',access='append')
		write(56,'(i10)') Nslice
	close(56)
	allocate(MolNumSl(SubNum+1,Nslice+1))
	allocate(Vol1MolNum(SubNum+1))
	allocate(Vol2MolNum(SubNum+1))
	allocate(MemIn1(SubNum+1))
	allocate(MemOut1(SubNum+1))
	allocate(InMem1(SubNum+1))
	allocate(MemIn2(SubNum+1))
	allocate(MemOut2(SubNum+1))
	allocate(InMem2(SubNum+1))
	
	if (First==1) then
		First=0
		do i=1,SubNum
		MemIn1(i)=0.0
		MemOut1(i)=0.0
		InMem1(i)=0.0
		MemIn2(i)=0.0
		MemOut2(i)=0.0
		InMem2(i)=0.0
		NearInVol(i)=0.0
		NearOutVol(i)=0.0
		Vol1MolNum=0.0
		Vol2MolNum=0.0
			do j=1,Nslice
				MolNumSl(i,j)=0.0
			enddo
		enddo
		SampleNum=0
		open(28,file='SlDens.out')
			write(28,'(a)')  'r,(nm)       ntbef,(nm-3)       ntaf,(nm-3)'
		close(28)
		open(27,file='TimeDep.out')
			write(27,'(a,$)') 'time,(ps)   '
			!do i=1,SubNum
			!	write(TempString,'(i5)') SubNum
			!	write(27,'(3a,$)') ' ak(',trim(adjustl(TempString)),'),x  '
			!enddo
			do i=1,SubNum
				write(27,'(3a,$)') '  MemIn1(',trim(adjustl(TempString)),') '
				write(27,'(3a,$)') '  InMem1(',trim(adjustl(TempString)),') '
				write(27,'(3a,$)') '  MemOut1(',trim(adjustl(TempString)),') '
				write(27,'(3a,$)') '  MemIn2(',trim(adjustl(TempString)),') '
				write(27,'(3a,$)') '  InMem2(',trim(adjustl(TempString)),') '
				write(27,'(3a,$)') '  MemOut2(',trim(adjustl(TempString)),') '
			enddo
			do i=1,SubNum
				write(27,'(3a,$)') '  Vol1Dens',trim(adjustl(TempString)),'(nm-1)  '
				write(27,'(3a,$)') '  Vol1Dens*',trim(adjustl(TempString)), ' '
				write(27,'(3a,$)') '  Vol2Dens',trim(adjustl(TempString)),'(nm-1)  '
				write(27,'(3a,$)') '  Vol2Dens*',trim(adjustl(TempString)), ' '
			enddo
			do i=1,SubNum
				write(27,'(3a,$)')  'Flow1',trim(adjustl(TempString)),'(nm-2ps-1)  '
				write(27,'(3a,$)')  'Flow2',trim(adjustl(TempString)),'(nm-2ps-1)  '
			enddo
			do i=1,SubNum
				write(27,'(3a,$)')  'NF1',trim(adjustl(TempString)),'  '
				write(27,'(3a,$)')  'NF2',trim(adjustl(TempString)),'  '
			enddo
			do i=1,SubNum
				write(27,'(3a,$)')  'Ntot',trim(adjustl(TempString)),'  '
			enddo
			write(27,'(a)') ' '
		close(27)
	else
		open(28,file='SlDens.out')
			write(28,'(a)')  'r,(nm)       ntbef,(nm-3)       ntaf,(nm-3)'
		close(28)
		open(29,file='calc.temp')
		
		do i=1,SubNum
			MemIn1(i)=0.0
			MemOut1(i)=0.0
			InMem1(i)=0.0
			MemIn2(i)=0.0
			MemOut2(i)=0.0
			InMem2(i)=0.0
			Vol1MolNum=0.0
			Vol2MolNum=0.0
			do j=1,Nslice
				read(29,'(f20.5)')MolNumSl(i,j)
			enddo
			read(29,'(f20.5)') NearInVol(i)
			read(29,'(f20.5)') NearOutVol(i)
		enddo
		read(29,*) SampleNum
		close(29)
	endif
	
	SampleNum=SampleNum+1
	do i=1,SubNum
		do j=1,NLiq(i)
			hist=ceiling(float(Nslice)*Vol2MolZ(i,j)/TotBoxLen)
			MolNumSl(i,hist)=MolNumSl(i,hist)+1.0
			if ((Vol1MolZ(i,j)>BoxL1).and.(Vol1MolZ(i,j)<BoxL2)) then
				Vol1MolNum(i)=Vol1MolNum(i)+1.0
			endif
			if ((Vol2MolZ(i,j)>BoxL1).and.(Vol2MolZ(i,j)<BoxL2)) then
				Vol2MolNum(i)=Vol2MolNum(i)+1.0
			endif
			!
			if ((Vol1MolZ(i,j)>Mem1E).and.(Vol1MolZ(i,j)<Mem2B)) then
				MemIn1(i)=MemIn1(i)+1.0
			elseif ((Vol1MolZ(i,j)>Mem2E)) then
				MemOut1(i)=MemOut1(i)+1.0
			else
				InMem1(i)=InMem1(i)+1.0
			endif
			!
			if ((Vol2MolZ(i,j)>Mem1E).and.(Vol2MolZ(i,j)<Mem2B)) then
				MemIn2(i)=MemIn2(i)+1.0
			elseif ((Vol2MolZ(i,j)>Mem2E)) then
				MemOut2(i)=MemOut2(i)+1.0
			else
				InMem2(i)=InMem2(i)+1.0
			endif
			!
			if ((Vol2MolZ(i,j)<MemDelta*Sigma/2.0).or.(Vol2MolZ(i,j)>TotBoxLen-MemDelta*Sigma/2.0)) then
				NearOutVol(i)=NearOutVol(i)+1.0
			elseif ((Vol2MolZ(i,j)>Mem2E-MemDelta*Sigma/2.0).and.(Vol2MolZ(i,j)<Mem2E+MemDelta*Sigma/2.0)) then
				NearOutVol(i)=NearOutVol(i)+1.0
			endif
			if ((Vol2MolZ(i,j)>Mem1E-MemDelta*Sigma/2.0).and.(Vol2MolZ(i,j)<Mem1E+MemDelta*Sigma/2.0)) then
				NearInVol(i)=NearInVol(i)+1.0
			elseif ((Vol2MolZ(i,j)>Mem2B-MemDelta*Sigma/2.0).and.(Vol2MolZ(i,j)<Mem2B+MemDelta*Sigma/2.0)) then
				NearInVol(I)=NearInVol(i)+1.0
			endif
		enddo
	enddo
	open(34,file='calc.temp')
		do i=1,SubNum
			do j=1,Nslice
				write(34,'(f20.5)') MolNumSl(i,j)
			enddo
			write(34,'(f20.5)') NearInVol(i)
			write(34,'(f20.5)') NearOutVol(i)
		enddo
	write(34,'(i10)') SampleNum
	close(34)
	
	open(28,file='SlDens.out',access='append')
		do i=1,Nslice
			write(28,'(f10.5,$)') TotBoxLen*float(i)/float(Nslice)
			do j=1,SubNum
				write(28,'(2f30.5,$)') MolNumSl(j,i)/(BoxH*BoxH*TotBoxLen/float(Nslice))/float(SampleNum),&
				&MolNumSl(j,i)/(BoxH*BoxH*TotBoxLen/float(Nslice))*Sigma*Sigma*Sigma/float(SampleNum)
			enddo
			write(28,'(a)') ' '
		enddo
	close(28)
	
	open(39,file='NearDens.out')
	do i=1,SubNum
		write(39,'(5f30.5)') NearInVol(i)/(BoxH*BoxH*MemDelta*Sigma)/float(SampleNum)/2.0,&
		&NearInVol(i)/(BoxH*BoxH*MemDelta*Sigma)/float(SampleNum)*Sigma*Sigma*Sigma/2.0,&
		&NearOutVol(i)/(BoxH*BoxH*MemDelta*Sigma)/float(SampleNum)/2.0,&
		&NearOutVol(i)/(BoxH*BoxH*MemDelta*Sigma)/float(SampleNum)*Sigma*Sigma*Sigma/2.0,&
		&(NearInVol(i)-NearOutVol(i))/(BoxH*BoxH*MemDelta*Sigma)/float(SampleNum)/(Mem1E-Mem1B+MemDelta*Sigma)/2.0
	enddo
	close(39)
	
	open(29,file='TimeDep.out',access='append')
		write(29,'(f20.10,$)') CurTime-0.5*DTime
!		do i=1,SubNum
!			write(29,'(f20.10,$)') AkL(i)
!		enddo
		do i=1,SubNum
			write(29,'(6f20.10,$)') MemIn1(i),InMem1(i),MemOut1(i),MemIn2(i),InMem1(i),MemOut2(i)
		enddo
		do i=1,SubNum
			write(29,'(f20.10,$)') Vol1MolNum(i)/(BoxH*BoxH*BoxH)
			write(29,'(f20.10,$)') Vol1MolNum(i)/(BoxH*BoxH*BoxH)*Sigma*Sigma*Sigma
			write(29,'(f20.10,$)') Vol2MolNum(i)/(BoxH*BoxH*BoxH)
			write(29,'(f20.10,$)') Vol2MolNum(i)/(BoxH*BoxH*BoxH)*Sigma*Sigma*Sigma
		enddo
		do i=1,SubNum
			write(29,'(f20.10,$)') (MemIn1(i)-MemIn2(i))/(BoxH*BoxH)/DTime
			write(29,'(f20.10,$)') (MemOut2(i)-MemOut1(i))/(BoxH*BoxH)/DTime 
		enddo
		do i=1,SubNum
			write(29,'(f20.10,$)') (MemIn1(i)-MemIn2(i))
			write(29,'(f20.10,$)') (MemOut2(i)-MemOut1(i))
		enddo
		do i=1,SubNum
			write(29,'(i7,$)') NLiq(i)
		enddo
		write(29,'(a)') ' '
	close(29)
	
	open(7,file='test.temp')
		write(7,'(f20.10)') BoxL1
		write(7,'(f20.10)') BoxL2
		write(7,'(f20.10)') BoxG1
		write(7,'(f20.10)') BoxG2
		write(7,'(f20.10)') BoxH
		write(7,'(i6)') SubNum
		write(7,'(i6)') TotalMol
		do i=1,SubNum
			write(7,'(a5,2i6)') SubName(i),NLiq(i),SubAtomNum(i)
		enddo
		write(7,'(i6)') MemAtom
		write(7,'(f20.10)') TotBoxLen
		write(7,'(f20.10)') Mem1B
		write(7,'(f20.10)') Mem1E
		write(7,'(f20.10)') Mem2B
		write(7,'(f20.10)') Mem2E
		write(7,'(f20.10)') CurTime
		write(7,'(i6)') First
	close(7)
	print *, 'MAIN IN FILE Write DONE'
	
end program

