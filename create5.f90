program mainv5
	implicit none
	integer(4) SubNum	!unmber of substances
	character(20) SubGroFile
	character(20) TempString
	integer(4),allocatable:: InitAtom(:)
	integer(4),allocatable:: InitNum(:,:)
	character(5),allocatable:: InitSubName(:,:)
	character(5),allocatable:: InitAtomName(:,:)
	real(8),allocatable:: InitX(:,:)
	real(8),allocatable:: InitY(:,:)
	real(8),allocatable:: InitZ(:,:)
	real(8),allocatable:: InitVX(:,:)
	real(8),allocatable:: InitVY(:,:)
	real(8),allocatable:: InitVZ(:,:)
	integer(4) MaxAtoms
	integer(4) i,j,k
	integer(4) TempInt
	integer(4) MemType
	integer(4) Mem1HW
	integer(4) Mem1Len
	real(8) Sigma
	real(8) RoL
	real(8),allocatable:: frac(:)
	integer(4),allocatable:: NLiq(:)
	real(8) BoxGasLen
	real(8) BoxH
	real(8) BoxW
	real(8) BoxLiqLen
	real(8) BoxLiqVol
	real(8) MemLen
	integer(4) NLiqTot
	integer(4) BoxPart
	real(8) MemDelta
	real(8),allocatable:: MemX(:),MemY(:),MemZ(:)
	real(8),allocatable:: BoxArrX(:),BoxArrY(:),BoxArrZ(:)
	integer(4), allocatable:: BoxArrOpen(:)
	integer(4),allocatable:: MolTypeNum(:)
	integer(4) RandInt
	integer(4) Pass
	integer(4) NAtomTot
	real(8) eps1
	real(8),allocatable:: Ak1(:),Ak2(:)
	integer(4) Nstep
	real(8) Temp
	real(8) CurTime,DTime
	integer(4) First
	!
	MaxAtoms=40
	call srand(15872)
	eps1=0.00831446
	
	!initial files 
	open(7,file='main.in')
		read(7,'(a)') TempString
		read(7,*) SubNum
		print *, TempString,SubNum
		allocate(InitAtom(SubNum))
		allocate(InitNum(SubNum,MaxAtoms))
		allocate(InitSubName(SubNum,MaxAtoms))
		allocate(InitAtomName(SubNum,MaxAtoms))
		allocate(InitX(SubNum,MaxAtoms))
		allocate(InitY(SubNum,MaxAtoms))
		allocate(InitZ(SubNum,MaxAtoms))
		allocate(InitVX(SubNum,MaxAtoms))
		allocate(InitVY(SubNum,MaxAtoms))
		allocate(initVZ(SubNum,MaxAtoms))
		
		allocate(frac(SubNum))
		allocate(NLiq(SubNum))
		allocate(Ak1(SubNum))
		allocate(Ak2(SubNum))
		
		do i=1,SubNum	!reading grofile
			read(7,*) SubGroFile, frac(i), Ak1(i),Ak2(i)
			write(SubGroFile,'(2a)') trim(adjustl(SubGroFile)), '.gro'
			open(8,file=SubGroFile)
			read(8, '(a)') TempString
			read(8,'(i6)') InitAtom(i)
			do j=1,InitAtom(i)
				read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, InitSubName(i,j), InitAtomName(i,j),&
				&InitNum(i,j),&
				&InitX(i,j),InitY(i,j),InitZ(i,j),&
				&InitVX(i,j),InitVY(i,j),InitVZ(i,j)
			enddo
			close(8)
			print *,' Reading ', SubGroFile, '  DONE'
			print *,' Fraction ', frac(i), ' Activiti L ', Ak1(i) ,' Activiti G ', Ak2(i) 
		enddo	!reading done
		read(7,'(a)') TempString
		read(7,*) MemType
		print*, TempString, MemType 
		if (MemType==1) then
			read(7,'(a)') TempString
			read(7,*) Sigma
			print *, TempString, Sigma
			read(7,'(a)') TempString
			read(7,*) MemDelta
			print*,TempString, MemDelta, ' sigmas'
			read(7,'(a)') TempString
			read(7,*) Mem1HW
			print *, TempString, Mem1HW
			read(7,'(a)') TempString
			read(7,*) Mem1Len
			print *, TempString, Mem1Len
			read(7,'(a)') TempString
			read(7,*) RoL
			print*,TempString,RoL
			read(7,'(a)') TempString
			read(7,*) Nstep
			print*,TempString,Nstep
			read(7,'(a)') TempString
			read(7,*) Temp
			print*,TempString,Temp
			read(7,'(a)') TempString
			read(7,*) Dtime
			if (SubNum>1) then
				
			endif
		endif
		close(7)
		if (MemType==2) then
			
		endif
		
		print*, 'main.in file DONE'
		print*,'  '
		
		if (MemType==1) then
			BoxH=Mem1HW*Sigma*MemDelta
			BoxW=Mem1HW*Sigma*MemDelta
			MemLen=Mem1Len*Sigma*MemDelta
			BoxLiqLen=3.0*BoxH
			BoxGasLen=3.0*BoxLiqLen
			BoxLiqVol=BoxH*BoxW*BoxLiqLen
			NLiqTot=RoL*BoxLiqVol/Sigma/Sigma/Sigma
			if (SubNum==1) then
				NLiq(1)=NLiqTot
			else
				TempInt=0
				do i=1,SubNum-1
					NLiq(i)=ceiling(NLiqTot*frac(i))
					TempInt=TempInt+Nliq(i)
				enddo
				Nliq(SubNum)=NLiqTot-TempInt
			endif
			NAtomTot=0
			do i=1,SubNum
				NAtomTot=NAtomTot+InitAtom(i)*Nliq(i)
			enddo
			BoxPart=1
			do while (BoxPart*BoxPart*BoxPart*3<NLiqTot)
				BoxPart=BoxPart+1
			enddo
			print *, ' NLiqTot ', NLiqTot
			allocate(MemX(Mem1HW*Mem1HW*Mem1Len))
			allocate(MemY(Mem1HW*Mem1HW*Mem1Len))
			allocate(MemZ(Mem1HW*Mem1HW*Mem1Len))
			allocate(BoxArrX(BoxPart*BoxPart*BoxPart*3))
			allocate(BoxArrY(BoxPart*BoxPart*BoxPart*3))
			allocate(BoxArrZ(BoxPart*BoxPart*BoxPart*3))
			allocate(BoxArrOpen(BoxPart*BoxPart*BoxPart*3))
			do i=1,Mem1Len
				do j=1,Mem1HW
					do k=1,Mem1HW
						MemX((i-1)*Mem1HW*Mem1HW+(j-1)*Mem1Hw+k)=(float(j)-0.5)*Sigma*MemDelta
						MemY((i-1)*Mem1HW*Mem1HW+(j-1)*Mem1Hw+k)=(float(k)-0.5)*Sigma*MemDelta
						MemZ((i-1)*Mem1HW*Mem1HW+(j-1)*Mem1Hw+k)=(float(i)-0.5)*Sigma*MemDelta
					enddo
				enddo
			enddo
			print *, ' Mem Creation DONE'
			do i=1,BoxPart*3
				do j=1,BoxPart
					do k=1,BoxPart
						BoxArrX((i-1)*BoxPart*BoxPart+(j-1)*BoxPart+k)=(float(j)-0.5)*BoxH/float(BoxPart)
						BoxArrY((i-1)*BoxPart*BoxPart+(j-1)*BoxPart+k)=(float(k)-0.5)*BoxW/float(BoxPart)
						BoxArrZ((i-1)*BoxPart*BoxPart+(j-1)*BoxPart+k)=(float(i)-0.5)*BoxLiqLen/float(BoxPart*3)
						BoxArrOpen((i-1)*BoxPart*BoxPart+(j-1)*BoxPart+k)=0
					enddo
				enddo
			enddo
			print *, ' BoxArr Creation DONE'
			allocate(MolTypeNum(NLiqTot))
			TempInt=0
			do i=1,SubNum
				do j=1,NLiq(i)
					Pass=0
					TempInt=TempInt+1
					do while (Pass==0)
						RandInt=ceiling(rand()*float(BoxPart*BoxPart*BoxPart*3))
						if (BoxArrOpen(RandInt)==0) then
							MolTypeNum(TempInt)=RandInt
							BoxArrOpen(RandInt)=i
							pass=1
						endif
					enddo
				enddo
			enddo
			print *, ' Randomization BoxArr DONE '
		endif
		!generating GRO
		print *, NLiqTot
		open(20,file='test.gro')
			write(20,'(a)') ' generated by memcreat v5'
			write(20,'(i6)') NAtomTot+Mem1HW*Mem1HW*Mem1Len*2 !BoxPart*BoxPart*BoxPart*3+
			TempInt=1
			do i=1,NLiqTot
				do j=1,InitAtom(BoxArrOpen(MolTypeNum(i)))
					write(20,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, InitSubName(BoxArrOpen(MolTypeNum(i)),j),&
					&InitAtomName(BoxArrOpen(MolTypeNum(i)),j),&
					&InitNum(BoxArrOpen(MolTypeNum(i)),j),&
					&InitX(BoxArrOpen(MolTypeNum(i)),j)+BoxArrX(MolTypeNum(i)),&
					&InitY(BoxArrOpen(MolTypeNum(i)),j)+BoxArrY(MolTypeNum(i)),&
					&InitZ(BoxArrOpen(MolTypeNum(i)),j)+BoxArrZ(MolTypeNum(i))+MemLen,&
					&InitVX(BoxArrOpen(MolTypeNum(i)),j),&
					&InitVY(BoxArrOpen(MolTypeNum(i)),j),&
					&InitVZ(BoxArrOpen(MolTypeNum(i)),j)
					TempInt=TempInt+1
				enddo
			enddo

!			do i=1,BoxPart*BoxPart*BoxPart*3
!				write(20,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, ' B ',&
!					&'  B ',&
!					& 1 ,&
!					&BoxArrX(i),&
!					&BoxArrY(i),&
!					&BoxArrZ(i)+MemLen,&
!					&0.0,0.0,0.0
!					TempInt=TempInt+1
!			enddo
			do i=1,Mem1HW*Mem1HW*Mem1Len
				write(20,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, 'MEM', ' A',&
				&TempInt,MemX(i), MemY(i),MemZ(i),&
				&0.0,0.0,0.0
				TempInt=TempInt+1
			enddo
			do i=1,Mem1HW*Mem1HW*Mem1Len
				write(20,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, 'MEM', ' A',&
				&TempInt,&
				&MemX(i),MemY(i),MemZ(i)+MemLen+BoxLiqLen,&
				&0.0,0.0,0.0
				TempInt=TempInt+1
			enddo
			write(20,'(3f20.5)') BoxH, BoxW, MemLen*2+BoxLiqLen+BoxGasLen
		close(20)
		!top files
		open(21,file='test.top')
			write(21,'(a)') ' [ defaults ] '
			write(21,'(2i5,a,2f10.5)') 1, 2, ' no ', 1.0 , 1.0

			write(21,'(a)') ' [ atomtypes ] '
			write(21,'(a,2f10.4,a,2e20.10)') ' O ', 1.0, 0.0, ' A ', Sigma*1.0,  eps1 
			write(21,'(a,2f10.4,a,2e20.10)') ' A ', 1.0, 0.0, ' A ', Sigma*1.0,  eps1 
			write(21,'(a,2f10.4,a,2e20.10)') ' B ', 1.0, 0.0, ' A ', Sigma*1.0,  eps1 
			write(21,'(a,2f10.4,a,2e20.10)') ' C ', 1.0, 0.0, ' A ', Sigma*1.0,  eps1 
			
			write(21,'(a)') ' [ moleculetype ] '
			write(21,'(a,a)') '  B   ',  ' 0'
			write(21,'(a)') ' [ atoms ] '
			write(21,'(i5,a,i5,2a,i5,2f10.5)') 1 ,' B ',1,' B ',' B ',1,0.00,1.0 
			write(21,'(a)') ' [ bonds ] '
			write(21,'(a)') ' [ pairs ] '
			write(21,'(a)') ' [ angles ] '
			write(21,'(a)') ' [ dihedrals ] '
			write(21,'(a)') ' [ exclusions ] '
			write(21,'(a)') ' [ position_restraints ] '
			write(21,'(a)') ' [ constraints ] '
			write(21,'(a)') '    '
			
			write(21,'(a)') ' [ moleculetype ] '
			write(21,'(a,a)') '  C   ',  ' 0'
			write(21,'(a)') ' [ atoms ] '
			write(21,'(i5,a,i5,2a,i5,2f10.5)') 1 ,' C ',1,' C ',' C ',1,0.00,1.0 
			write(21,'(a)') ' [ bonds ] '
			write(21,'(a)') ' [ pairs ] '
			write(21,'(a)') ' [ angles ] '
			write(21,'(a)') ' [ dihedrals ] '
			write(21,'(a)') ' [ exclusions ] '
			write(21,'(a)') ' [ position_restraints ] '
			write(21,'(a)') ' [ constraints ] '
			write(21,'(a)') '    '
			
			write(21,'(a)') ' [ moleculetype ] '
			write(21,'(a,a)') '  MEM ',  ' 0'
			write(21,'(a)') ' [ atoms ] '
			write(21,'(i5,a,i5,2a,i5,2f10.5)') 1 ,' A ',1,' A ',' A ',1,0.00,1.0 
			write(21,'(a)') ' [ bonds ] '
			write(21,'(a)') ' [ pairs ] '
			write(21,'(a)') ' [ angles ] '
			write(21,'(a)') ' [ dihedrals ] '
			write(21,'(a)') ' [ exclusions ] '
			write(21,'(a)') ' [ position_restraints ] '
			write(21,'(a)') ' [ constraints ] '
			write(21,'(a)') '    '
			
			write(21,'(a)') ' [ system ] '
			write(21,'(a)') ' generateg with membcreat v5 '
			write(21,'(a)') ' [ molecules ] '
			do i=1,SubNum
				write(21,'(a,i6)') InitSubName(i,1), NLiq(i)
			enddo
			write(21,'(a,i6)') ' MEM  ', Mem1HW*Mem1HW*Mem1Len*2
		close(21)
		
		open(22,file='index.ndx')
			write(22,'(a)') '[ System ]'
			TempInt=1
			do i=1,NAtomTot+Mem1HW*Mem1HW*Mem1Len*2
				write(22,'(i6,$)') TempInt
				if (mod(TempInt,10)==0) then
					write(22,'(a)') '  '
				endif
				TempInt=TempInt+1
			enddo
			if (mod(TempInt,10)/=0) then
				write(22, '(a)') '  '
			endif
			write(22,'(a)') '[ B ] ' 
			TempInt=1
			do i=1,SubNum
				do j=1,Nliq(i)
					do k=1,InitAtom(i)
						write(22,'(i6,$)') TempInt
						if (mod(TempInt,10)==0) then
							write(22,'(a)') '  '
						endif
						TempInt=TempInt+1
					enddo
				enddo
			enddo
			if (mod(TempInt,10)/=0) then
				write(22, '(a)') '  '
			endif
			write(22, '(a)') '[ MEM ]'
			do i=1,Mem1HW*Mem1HW*Mem1Len*2
				write(22,'(i5,$)') TempInt
				if (mod(TempInt,10)==0) then
					write(22,'(a)') '  '
				endif
				TempInt=TempInt+1
			enddo
			if (mod(TempInt,10)/=0) then
				write(22, '(a)') '  '
			endif
		close(22)
		
		open(23,file='start.sh')
			write(23,'(a)') '#/bin/bash'
			write(23,'(a)') 'rm test.trp'
			write(23,'(a)') 'for n in {1..1}'
			write(23,'(a)') 'do'
				write(23,'(a)') 'echo $n " step" '
				write(23,'(a)') './dcvmd'
				write(23,'(a)') 'rm test.gro'
				write(23,'(a)') 'rm test.top'
				write(23,'(a)') 'rm index.ndx'
				write(23,'(a)') 'cp testout.gro test.gro'
				write(23,'(a)') 'cp testnew.top test.top'
				write(23,'(a)') 'cp indexnew.ndx index.ndx'
				write(23,'(a)') 'rm test.tpr'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') 'rm test.edr'
				write(23,'(a)') 'rm test.log'
				write(23,'(a)') 'rm test.trr'
				write(23,'(a)') 'rm test.xtc'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') "find -type f -name 'test.mdp' -print0 | xargs --null perl -pi -e &
				&'s/integrator               = md/integrator               = steep/'"
				write(23,'(a)') 'mygrompp -f test.mdp -c test.gro -n index.ndx -p test.top -o test'
				write(23,'(a)') 'mymdrun -deffnm test -pd'
				write(23,'(a)') 'cp test.gro test1.gro'
				write(23,'(a)') 'rm test.tpr'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') 'rm test.edr'
				write(23,'(a)') 'rm test.log'
				write(23,'(a)') 'rm test.trr'
				write(23,'(a)') 'rm test.xtc'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') './tempprog'
				write(23,'(a)') 'cp testnew.gro test.gro '
				write(23,'(a)') "find -type f -name 'test.mdp' -print0 | xargs --null perl -pi -e &
				&'s/integrator               = steep/integrator               = md/'"
				write(23,'(a)') 'mygrompp -f test.mdp -c test.gro -n index.ndx -p test.top -o test '
				write(23,'(a)') 'mymdrun -deffnm test -pd'
				write(23,'(a)') 'cp test.gro test2.gro '
				write(23,'(a)') 'rm test.tpr'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') 'rm test.edr'
				write(23,'(a)') 'rm test.log'
				write(23,'(a)') 'rm test.trr'
				write(23,'(a)') 'rm test.xtc'
				write(23,'(a)') 'rm test.cpt'
				write(23,'(a)') './denscalc'
				write(23,'(a)') 'rm \#*'
!				write(23,)
			write(23,'(a)') 'done'
		close(23)
		CurTime=0.0
		First=1
		open(24,file='test.temp')
			write(24,'(f20.10)') MemLen+BoxH	!liq vol 1
			write(24,'(f20.10)') MemLen+2.0*BoxH	!liq vol 2
			write(24,'(f20.10)') 2.0*MemLen+BoxLiqLen*2.0	!gas vol 1
			write(24,'(f20.10)') 2.0*MemLen+BoxLiqLen*3.0	!gas vol 2
			write(24,'(f20.10)') BoxH
			write(24,'(i6)') SubNum
			write(24,'(i6)') NLiqTot
			do i=1,SubNum
				write(24,'(a5,i6,i6)') InitSubName(i,1), NLiq(i), InitAtom(i)
			enddo
			write(24,'(i6)') Mem1HW*Mem1HW*Mem1Len*2
			write(24,'(f20.10)') MemLen*2+BoxLiqLen+BoxGasLen
			write(24,'(f20.10)') 0.0	!мембрана 1 начало
			write(24,'(f20.10)') MemLen	!мембрана 1 конец
			write(24,'(f20.10)') MemLen+BoxLiqLen !мембрана 2 начало
			write(24,'(f20.10)') MemLen*2+BoxLiqLen	!мемебрана 2 конец
			write(24,'(f20.10)') CurTime
			write(24,'(i5)') First
		close(24)

end program

