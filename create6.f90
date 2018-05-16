module GlobalVar
	implicit none
	integer(4) MemN
	character(20) SysName
	character(20) TempString
	integer(4) TempInt
	integer(4),allocatable:: MemMol(:)
	character(3),allocatable:: MemName(:)
	real(8),allocatable:: MemX(:)	!coordinates
	real(8),allocatable:: MemY(:)
	real(8),allocatable:: MemZ(:)
	real(8),allocatable:: MemVX(:)	!velosities
	real(8),allocatable:: MemVY(:)
	real(8),allocatable:: MemVZ(:)
	real(8) MemHX	!box lenght
	real(8) MemHY
	real(8) MemHZ

	real(8) ConVol1B	!control voloume 1 begin
	real(8) ConVol1E	!end
	real(8) ConVol2B
	real(8) ConVol2E

	integer(4) SubNum
	integer(4),allocatable:: InitNum(:,:)
	integer(4),allocatable:: InitAtom(:)
	character(5),allocatable:: InitSubName(:,:)
	character(3),allocatable:: InitAtomName(:,:)
	real(8),allocatable:: InitX(:,:)
	real(8),allocatable:: InitY(:,:)
	real(8),allocatable:: InitZ(:,:)
	real(8),allocatable:: InitVX(:,:)
	real(8),allocatable:: InitVY(:,:)
	real(8),allocatable::InitVZ(:,:)
	integer(4) MaxAtoms
	character(20) SubGroFile
end module

program generate
	use GlobalVar
	implicit none

	MaxAtoms=20
	print *, 'Program init'
	call ReadMem()
	ConVol1B=MemHZ+MemHX
	ConVol1E=ConVol1B+MemHX
	ConVol2B=MemHZ*2.0+MemHX*3.0+MemHX*3.0
	ConVol2E=ConVol2B+MemHX*3.0

	call ReadMainIn()
	call WriteMem()
	call WriteTop()
	call WriteScript()
	call WriteTemp()

	print *, 'Program DONE'
end program

subroutine ReadMem( )
	use GlobalVar
	implicit none
	integer(4) i

	open(11,file='mem.gro')
		read(11,*) SysName
		read(11,*) MemN
		allocate(MemMol(MemN))
		allocate(MemName(MemN))
		allocate(MemX(MemN))
		allocate(MemY(MemN))
		allocate(MemZ(MemN))
		allocate(MemVX(MemN))
		allocate(MemVY(MemN))
		allocate(MemVZ(MemN))
		do i=1,MemN
			read(11,'(i5,2a5,i5,3f8.3,3f8.4)') MemMol(i), TempString, MemName(i),&
			&TempInt,&
			&MemX(i),MemY(i),MemZ(i),&
			&MemVX(i),MemVY(i),MemVZ(i)
		enddo
		read(11,*) MemHX, MemHY, MemHZ
	close(11)
	print *, 'membrane reading DONE'
end subroutine

subroutine WriteMem()
	use GlobalVar
	implicit none
	integer(4) i,j
	integer(4) SumAtomsOf2
	integer(4) tempi

	SumAtomsOf2=0
	do i=1,SubNum
		do j=1,InitAtom(i)
			SumAtomsOf2=SumAtomsOf2+1
		enddo
	enddo
	print *, SubNum,SumAtomsOf2
	tempi=1
	open(21,file='test.gro')
		write(21,'(a)') SysName
		write(21,'(i5)') MemN*2+SumAtomsOf2
		do i=1,SubNum
			do j=1,InitAtom(i)
				write(21,'(i5,2a5,i5,3f8.3,3f8.4)') i, InitSubName(i,j), InitAtomName(i,j),&
				&tempi,&
				&InitX(i,j),InitY(i,j),InitZ(i,j)+MemHZ+MemHX*(1.0+float(i)/10.0),&
				&0.0,0.0,0.0
				tempi=tempi+1
			enddo
		enddo
		do i=1,2
			do j=1,MemN
				write(21,'(i5,2a5,i5,3f8.3,3f8.4)') MemMol(j)+2+MemMol(MemN)*(i-1), 'MEM', MemName(j),&
			&tempi,&
			&MemX(j),MemY(j),MemZ(j)+float(i-1)*(MemHZ+3.0*MemHX),&
			&0.0,0.0,0.0
			tempi=tempi+1
			enddo
		enddo
		write(21,*) MemHX, MemHX, MemHZ*2.0+MemHX*(3.0+9.0)
	close(21)

	print *, 'Write gro file DONE'
end subroutine

subroutine WriteTop()
	use GlobalVar
	implicit none
	integer(4) i,j
	real(8) ChSi,ChC2,ChO
	real(8) MassSi,MassC2,MassO

	ChSi=0.9
	ChC2=-0.225
	ChO=-0.45

	MassSi=28.0
	MassC2=14.0
	MassO=16

	open(21,file='test.top')
		write(21,'(a)') ' [ defaults ] '
		write(21,'(2i5,a,2f10.5)') 1, 2, ' no ', 1.0 , 1.0

		write(21,'(a)') ' [ atomtypes ] '
		write(21,'(a,2f10.4,a,2f10.5)') ' Si ', MassSi, ChSi, ' A ', 0.38264,1.257
		write(21,'(a,2f10.4,a,2f10.5)') ' C2',  MassC2, ChC2, ' A ', 0.395,0.382444
		write(21,'(a,2f10.4,a,2f10.5)') ' O ',  MassO,  ChO,  ' A ', 0.311814,0.6285
		!ethanol
		write(21,'(a,2f10.4,a,2f10.5)') ' CE1 ',  15.011 , 0.000,   ' A ', 0.36072, 0.99898
		write(21,'(a,2f10.4,a,2f10.5)') ' CE2 ',  14.011 , 0.2526,  ' A ', 0.34612, 0.717461
		write(21,'(a,2f10.4,a,2f10.5)') ' OE1 ',  15.9994,-0.69711, ' A ', 0.31496, 0.707168
		write(21,'(a,2f10.4,a,2f10.5)') ' HE1 ',   1.0080, 0.44151, ' A ', 0.01,   0.0
		!water
		write(21,'(a,2f10.4,a,2f10.5)') ' OW ',  15.9994, 0.00000, ' A ', 0.31589, 0.774903
		write(21,'(a,2f10.4,a,2f10.5)') ' HW ',   1.0080, 0.5564,  ' A ', 0.01,   0.0
		write(21,'(a,2f10.4,a,2f10.5)') ' MW ',   0.0000, -1.1128, ' D ', 0.01,   0.0

		write(21,'(a)') ' [ bondtypes ] '
		write(21,'(2a,i5,2f15.5)') '  CE1  ','  CE2  ', 1 ,  0.19842,  500000.0
		write(21,'(2a,i5,2f15.5)') '  CE2  ','  OE1  ', 1 ,  0.171581,  500000.0
		write(21,'(2a,i5,2f15.5)') '  OE1  ','  HE1  ', 1 ,  0.095053,  500000.0
		write(21,'(2a,i5,2f15.5)') '  OW   ','  HW  ', 1 ,  0.09572,  500000.0
		write(21,'(2a,i5,2f15.5)') '  OW   ','  MW  ', 1 ,  0.01546,  500000.0

		write(21,'(a)') ' [ angletypes ] '
		write(21,'(3a,i5,2f10.5)') '  CE1  ','  CE2  ', '  OE1  ', 1 , 90.95 ,  500.0
		write(21,'(3a,i5,2f10.5)') '  CE2  ','  OE1  ', '  HE1  ', 1 , 106.368 ,  500.0
		write(21,'(3a,i5,2f10.5)') '  HW   ','  OW   ', '  HW   ', 1 , 104.52 ,  500.0
		write(21,'(3a,i5,2f10.5)') '  MW  ','  OW  ', '  HW   ', 1 , 52.26 ,  500.0

		write(21,'(a)') ' [ dihedraltypes ] '
		write(21,'(4a,i5,4f10.5)') '  CE1  ','  CE2  ', '  OE1  ','  HE1  ', 5 ,5.0,0.0,0.0,0.0

		write(21,'(a)') ' [ moleculetype ] '
		write(21,'(a,a)') '  MEM ',  ' 0'
		write(21,'(a)') ' [ atoms ] '
		do i=1,MemN
			write(21,'(i5,a,i5,2a,i5,$)') i ,MemName(i),1,' MEM ', MemName(i), mod(i-1,20)+1
			if(trim(adjustl(MemName(i)))=='Si') then
				write(21,'(2f10.5)') ChSi,MassSi
			endif
			if(trim(adjustl(MemName(i)))=='C2') then
				write(21,'(2f10.5)') ChC2,MassC2
			endif
			if(trim(adjustl(MemName(i)))=='O') then
				write(21,'(2f10.5)') ChO,MassO
			endif
		enddo
		write(21,'(a)') ' [ bonds ] '
		write(21,'(a)') ' [ pairs ] '
		write(21,'(a)') ' [ angles ] '
		write(21,'(a)') ' [ dihedrals ] '
		write(21,'(a)') ' [ exclusions ] '
		write(21,'(a)') ' [ position_restraints ] '
		write(21,'(a)') ' [ constraints ] '
		write(21,'(a)') '    '

		write(21,'(a)') ' [ moleculetype ] '
		write(21,'(a,a)') '  ETH ',  ' 3'
		write(21,'(a)') ' [ atoms ] '

		write(21,'(i5,a,i5,2a,i5,2f10.5)') 1 ,'  CE1  ',1,' ETH ', 'CE1' , 1 , 0.000 , 15.011
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 2 ,'  CE2  ',1,' ETH ', 'CE2' , 1 , 0.2556 , 14.011
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 3 ,'  OE1  ',1,' ETH ', 'OE1' , 1 , -0.69711 , 15.9994
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 4 ,'  HE1  ',1,' ETH ', 'HE1' , 1 , 0.44151 , 1.008

		write(21,'(a)') ' [ bonds ] '

		write(21,'(a)') ' 1    2    1'
		write(21,'(a)') ' 2    3    1'
		write(21,'(a)') ' 3    4    1'

		write(21,'(a)') ' [ pairs ] '
		write(21,'(a)') ' [ angles ] '

		write(21,'(a)') ' 1    2    3    1'
		write(21,'(a)') ' 2    3    4    1'

		write(21,'(a)') ' [ dihedrals ] '
		write(21,'(a)') '1    2    3    4    5'
		write(21,'(a)') ' [ exclusions ] '
		write(21,'(a)') ' [ position_restraints ] '
		write(21,'(a)') ' [ constraints ] '
		write(21,'(a)') '    '

		write(21,'(a)') ' [ moleculetype ] '
		write(21,'(a,a)') '  WAT ',  ' 3'
		write(21,'(a)') ' [ atoms ] '

		write(21,'(i5,a,i5,2a,i5,2f10.5)') 1 ,' OW ',1,' WAT ', ' OW ' , 1 , 0.000 , 15.9994
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 2 ,' HW ',1,' WAT ', ' HW ' , 1 , 0.5564 , 1.008  ! 1.0080, 0.5564
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 3 ,' HW ',1,' WAT ', ' HW ' , 1 , 0.5564 , 1.008
		write(21,'(i5,a,i5,2a,i5,2f10.5)') 4 ,' MW ',1,' WAT ', ' MW ' , 1 , -1.1128 , 0.0

		write(21,'(a)') ' [ bonds ] '

		write(21,'(a)') ' 1    2    1'
		write(21,'(a)') ' 1    3    1'

		write(21,'(a)') ' [ pairs ] '
		write(21,'(a)') ' [ angles ] '

		write(21,'(a)') ' 2    1    3    1'

		write(21,'(a)') ' [ dihedrals ] '
		write(21,'(a)') ' [ exclusions ] '
		write(21,'(a)') ' 1    2    3    4'
		write(21,'(a)') ' 2    1    3    4'
		write(21,'(a)') ' 3    1    2    4'
		write(21,'(a)') ' 4    1    2    3'

		!write(21,'(a)') ' [ position_restraints ] '
		!write(21,'(a)') ' [ constraints ] '
		write(21,'(a)') ' [ virtual_sites3 ] '
		write(21,'(a)') ' 4    1    2    3    1    0.131937768    0.131937768'
		write(21,'(a)') '    '
		write(21,'(a)') ' [ system ] '
		write(21,'(a)') ' generateg with membcreat v5 '
		write(21,'(a)') ' [ molecules ] '
		do i=1,SubNum
			write(21,'(a,i6)') InitSubName(i,1),1
		enddo
		write(21,'(a,i6)') ' MEM  ', 2

	close(21)
end subroutine

subroutine WriteScript()
	use GlobalVar
	implicit none
	integer(4) i

	open(23,file='start.sh')
		write(23,'(a)') '#/bin/bash'
		write(23,'(a)') 'rm test.trp'
		write(23,'(a)') 'for n in {1..1000}'
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
			write(23,'(a)') 'grompp -f test.mdp -c test.gro -n index.ndx -p test.top -o test -maxwarn 1 '
			write(23,'(a)') 'mdrun -deffnm test -pd'
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
			write(23,'(a)') 'grompp -f test.mdp -c test.gro -n index.ndx -p test.top -o test -maxwarn 1 '
			write(23,'(a)') 'mdrun -deffnm test -pd'
			write(23,'(a)') 'cp test.gro test2.gro '
			write(23,'(a)') 'rm test.tpr'
			write(23,'(a)') 'rm test.cpt'
			write(23,'(a)') 'rm test.edr'
			write(23,'(a)') 'rm test.log'
			write(23,'(a)') 'rm test.trr'
			write(23,'(a)') 'rm test.xtc'
			write(23,'(a)') 'rm test.cpt'
			write(23,'(a)') './dcvcalc'
			write(23,'(a)') 'rm \#*'
		write(23,'(a)') 'done'
	print *, ' Script creaton file DONE'

end subroutine

subroutine WriteTemp()
	use GlobalVar
	implicit none
	integer(4) i

	open(24,file='test.temp')
		write(24,'(f20.10)') ConVol1B	!liq vol 1
		write(24,'(f20.10)') ConVol1E	!liq vol 2
		write(24,'(f20.10)') ConVol2B	!gas vol 1
		write(24,'(f20.10)') ConVol2E	!gas vol 2
		write(24,'(f20.10)') MemHX
		write(24,'(i6)') SubNum	!number of substaces
		write(24,'(i6)') 0 !NLiqTot	!total liquid molecules/atoms ?
		do i=1,SubNum
			write(24,'(a5,i6,i6)') InitSubName(i,1),1,InitAtom(i) ! InitSubName(i,1), NLiq(i), InitAtom(i)
		enddo
		write(24,'(i6)') MemN*2	!Membrane atoms
		write(24,'(f20.10)') MemHZ*2.0+MemHX*(3.0+9.0)	!total box lenght
		write(24,'(f20.10)') 0.0	!мембрана 1 начало
		write(24,'(f20.10)') MemHZ	!мембрана 1 конец
		write(24,'(f20.10)') MemHZ+MemHX*3.0 !мембрана 2 начало
		write(24,'(f20.10)') MemHZ*2.0+MemHX*3.0	!мемебрана 2 конец
		write(24,'(f20.10)') 0.0 !CurTime	!current time
		write(24,'(i5)') 1 !First	!1 for the statistic to ero
	close(24)
	print *, ' Tempfile write DONE'
end subroutine

subroutine ReadMainIn()
	use GlobalVar
	implicit none
	integer(4) i,j
	real(8) TempReal1,TempReal2,TempReal3
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

		do i=1,SubNum	!reading grofile
			read(7,*) SubGroFile, TempReal1, TempReal2,TempReal3
			write(SubGroFile,'(2a)') trim(adjustl(SubGroFile)), '.gro'
			open(8,file=SubGroFile)
			read(8, '(a)') TempString
			read(8,'(i6)') InitAtom(i)
			do j=1,InitAtom(i)
				read(8,'(i5,2a5,i5,3f8.3,3f8.4)') TempInt, InitSubName(i,j), InitAtomName(i,j),&
				&InitNum(i,j),&
				&InitX(i,j),InitY(i,j),InitZ(i,j)
			enddo
			close(8)
			print *,' Reading ', SubGroFile, '  DONE'
			print *,' Fraction ', TempReal1, ' Activiti L ', TempReal2 ,' Activiti G ', TempReal3
		enddo	!reading done
		close(7)
end subroutine
