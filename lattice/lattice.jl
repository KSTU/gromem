"""
Generater gro and top files for gromacs calculation
-t --- type of latice [sc,bcc,fcc]
\t sc --- simple cubic
\t bcc --- body-centered cubic
\t fcc --- face-centered cubic 


julia lattice.jl -t sc -d 0.8 -l 0.125 -x 0.5
"""

using Printf


#get argumets
println("$(ARGS)")
id=0
for i in ARGS
    global id+=1
    if(i=="-t")
        global mType=ARGS[id+1]
    end
    if(i=="-d")
        global dens=parse(Float64,ARGS[id+1])
    end
    if(i=="-l")
        global latice_dens=parse(Float64,ARGS[id+1])
    end
    if(i=="-x")
        global xa=parse(Float64,ARGS[id+1])
    end
end

maxAtom=1500
sigma=0.3

tempV=maxAtom/dens  #box volume
tempL=tempV^(1.0/3.0)   #box length
   #numer of sites in membrane
#
sc_l=(1.0/latice_dens)^(1.0/3.0)
bcc_l=(2.0/latice_dens)^(1.0/3.0)
fcc_l=(4.0/latice_dens)^(1.0/3.0)


#membrane atoms




#cube
if(mType=="sc")
    n_size=(ceil(Int64,tempL/sc_l)+1)::Int64
    initText="simple cubic"
#    println(fileGro,"generated by latice.jl")
    #membrane atoms
    Nmem=n_size*n_size*n_size
    Len=n_size*sc_l
    Vol=Len^3
    NTotMol=ceil(Int64,dens*Vol)
    NAMol=ceil(Int64,xa*NTotMol)::Int64
    NBMol=(NTotMol-NAMol)::Int64
    #
    lbox_size=(1.0/(latice_dens*1.5))^(1.0/3.0)
    n2_size=ceil(Len/(1.0/(latice_dens*1.5))^(1.0/3.0))+1
    local n_cur::Int64=ceil(n2_size/n_size)*n_size+1
    println("n_cur $(n_cur) ncur^3 $(n_cur*n_cur*n_cur)")
    #get temp box
    temp_x=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_y=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_z=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_busy=Array{Float64}(undef,n_cur*n_cur*n_cur)
    print("num of places $(n_cur^3) num of molecules $(NTotMol)")
    #
    l_delta=Len/n_cur
    id=0
    for i=1:n_cur
        for j=1:n_cur
            for k=1:n_cur
                global id+=1
#                println("id $(id)  $(i*l_delta)")
                global temp_x[id]=i*l_delta
                global temp_y[id]=j*l_delta
                global temp_z[id]=k*l_delta
                global temp_busy[id]=0
            end
        end
    end
    println("length $(length(temp_x))")
    
    a_x=Array{Float64}(undef,NAMol)
    a_y=Array{Float64}(undef,NAMol)
    a_z=Array{Float64}(undef,NAMol)
    
    b_x=Array{Float64}(undef,NBMol)
    b_y=Array{Float64}(undef,NBMol)
    b_z=Array{Float64}(undef,NBMol)

    InsA=0
    InsB=0
    global id=0
    while(InsA<NAMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
#            println("rn $(rn) id  $(id)")
            global InsA+=1
            global id+=1
            a_x[id]=temp_x[rn]
            a_y[id]=temp_y[rn]
            a_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    global id=0
    while(InsB<NBMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
            global InsB+=1
            global id+=1
            b_x[id]=temp_x[rn]
            b_y[id]=temp_y[rn]
            b_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    
    id=0
    mem_x=Array{Float64}(undef,Nmem)
    mem_y=Array{Float64}(undef,Nmem)
    mem_z=Array{Float64}(undef,Nmem)
    for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-1)*sc_l
                mem_y[id]=(j-1)*sc_l
                mem_z[id]=(k-1)*sc_l
            end
        end
    end
end
#ock
if(mType=="bcc")
    n_size=(ceil(Int64,tempL/bcc_l)+1)::Int64
    initText="bcc"
#    println(fileGro,"generated by latice.jl")
    #membrane atoms
    Nmem=n_size*n_size*n_size*2	#2 latice
    Len=n_size*bcc_l
    Vol=Len^3
    NTotMol=ceil(Int64,dens*Vol)
    NAMol=ceil(Int64,xa*NTotMol)::Int64
    NBMol=(NTotMol-NAMol)::Int64
    #
    lbox_size=(1.0/(latice_dens*1.5))^(1.0/3.0)
    n2_size=ceil(Int64,Len/lbox_size)+1
#    n2_size=ceil(Int64,1.5*Len/bcc_l)+1
    n_cur=(ceil(Int64,n2_size/n_size*1.3)*n_size+1)::Int64
    println("n_cur $(n_cur) ncur^3 $(n_cur*n_cur*n_cur) $(n2_size) lboxsaiz $(lbox_size) Len $(Len) nsize $(n_size)")
    #get temp box
    temp_x=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_y=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_z=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_busy=Array{Float64}(undef,n_cur*n_cur*n_cur)
    print("num of places $(n_cur^3) num of molecules $(NTotMol)")
    #
    l_delta=Len/n_cur
    id=0
    for i=1:n_cur
        for j=1:n_cur
            for k=1:n_cur
                global id+=1
#                println("id $(id)  $(i*l_delta)")
                global temp_x[id]=i*l_delta
                global temp_y[id]=j*l_delta
                global temp_z[id]=k*l_delta
                global temp_busy[id]=0
            end
        end
    end
    println("length $(length(temp_x))")
    
    a_x=Array{Float64}(undef,NAMol)
    a_y=Array{Float64}(undef,NAMol)
    a_z=Array{Float64}(undef,NAMol)
    
    b_x=Array{Float64}(undef,NBMol)
    b_y=Array{Float64}(undef,NBMol)
    b_z=Array{Float64}(undef,NBMol)

    InsA=0
    InsB=0
    global id=0
    while(InsA<NAMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
#            println("rn $(rn) id  $(id)")
            global InsA+=1
            global id+=1
            a_x[id]=temp_x[rn]
            a_y[id]=temp_y[rn]
            a_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    global id=0
    while(InsB<NBMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
            global InsB+=1
            global id+=1
            b_x[id]=temp_x[rn]
            b_y[id]=temp_y[rn]
            b_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    
    id=0
    mem_x=Array{Float64}(undef,Nmem)
    mem_y=Array{Float64}(undef,Nmem)
    mem_z=Array{Float64}(undef,Nmem)
    for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-1)*bcc_l
                mem_y[id]=(j-1)*bcc_l
                mem_z[id]=(k-1)*bcc_l
            end
        end
    end
     for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-0.5)*bcc_l
                mem_y[id]=(j-0.5)*bcc_l
                mem_z[id]=(k-0.5)*bcc_l
            end
        end
    end
end

#gck
if(mType=="fcc")
    n_size=(ceil(Int64,tempL/fcc_l)+1)::Int64
    initText="fcc"
#    println(fileGro,"generated by latice.jl")
    #membrane atoms
    Nmem=n_size*n_size*n_size*4	#2 latice
    Len=n_size*fcc_l
    Vol=Len^3
    NTotMol=ceil(Int64,dens*Vol)
    NAMol=ceil(Int64,xa*NTotMol)::Int64
    NBMol=(NTotMol-NAMol)::Int64
    #
    lbox_size=(4.0/(latice_dens*1.5))^(1.0/3.0)
    n2_size=ceil(Len/(4.0/(latice_dens*1.5))^(1.0/3.0))+1
    local n_cur::Int64=ceil(n2_size/n_size*1.8)*n_size+1
    println("n_cur $(n_cur) ncur^3 $(n_cur*n_cur*n_cur)")
    #get temp box
    temp_x=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_y=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_z=Array{Float64}(undef,n_cur*n_cur*n_cur)
    temp_busy=Array{Float64}(undef,n_cur*n_cur*n_cur)
    print("num of places $(n_cur^3) num of molecules $(NTotMol)")
    #
    l_delta=Len/n_cur
    id=0
    for i=1:n_cur
        for j=1:n_cur
            for k=1:n_cur
                global id+=1
#                println("id $(id)  $(i*l_delta)")
                global temp_x[id]=i*l_delta
                global temp_y[id]=j*l_delta
                global temp_z[id]=k*l_delta
                global temp_busy[id]=0
            end
        end
    end
    println("length $(length(temp_x))")
    
    a_x=Array{Float64}(undef,NAMol)
    a_y=Array{Float64}(undef,NAMol)
    a_z=Array{Float64}(undef,NAMol)
    
    b_x=Array{Float64}(undef,NBMol)
    b_y=Array{Float64}(undef,NBMol)
    b_z=Array{Float64}(undef,NBMol)

    InsA=0
    InsB=0
    global id=0
    while(InsA<NAMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
#            println("rn $(rn) id  $(id)")
            global InsA+=1
            global id+=1
            a_x[id]=temp_x[rn]
            a_y[id]=temp_y[rn]
            a_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    global id=0
    while(InsB<NBMol)
        rn=rand(1:n_cur*n_cur*n_cur);
        if(temp_busy[rn]==0)
            global InsB+=1
            global id+=1
            b_x[id]=temp_x[rn]
            b_y[id]=temp_y[rn]
            b_z[id]=temp_z[rn]
            temp_busy[rn]=1
        end
    end
    
    id=0
    mem_x=Array{Float64}(undef,Nmem)
    mem_y=Array{Float64}(undef,Nmem)
    mem_z=Array{Float64}(undef,Nmem)
    for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-1)*fcc_l
                mem_y[id]=(j-1)*fcc_l
                mem_z[id]=(k-1)*fcc_l
            end
        end
    end
     for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-0.5)*fcc_l
                mem_y[id]=(j-0.5)*fcc_l
                mem_z[id]=(k-1)*fcc_l
            end
        end
    end
    for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-1)*fcc_l
                mem_y[id]=(j-0.5)*fcc_l
                mem_z[id]=(k-0.5)*fcc_l
            end
        end
    end
    for i=1:n_size
        for j=1:n_size
            for k=1:n_size
                global id+=1
                mem_x[id]=(i-0.5)*fcc_l
                mem_y[id]=(j-1)*fcc_l
                mem_z[id]=(k-0.5)*fcc_l
            end
        end
    end
end


#dimensionless
a_x=a_x*sigma
a_y=a_y*sigma
a_z=a_z*sigma

b_x=b_x*sigma
b_y=b_y*sigma
b_z=b_z*sigma

mem_x=mem_x*sigma
mem_y=mem_y*sigma
mem_z=mem_z*sigma
Len=Len*sigma

#save gro file
fileGro=open("test.gro","w")
println(fileGro,initText)
#println(fileGro,"temp string")
println(fileGro,"",NAMol+NBMol+Nmem)
for i=1:NAMol
    @printf(fileGro,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",i,"A","A",1,a_x[i],a_y[i],a_z[i])
end
for i=1:NBMol
    @printf(fileGro,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",i+NAMol,"B","B",1,b_x[i],b_y[i],b_z[i])
end
for i=1:Nmem
    @printf(fileGro,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",i+NAMol+NBMol,"M","M",1,mem_x[i],mem_y[i],mem_z[i])
end
#print(fileGro,"     $(Len)     $(Len)     $(Len)")
@printf(fileGro,"  %8.4f    %8.4f    %8.4f",Len,Len,Len)


#println("n_size $(n_size) sc_l $(sc_l) tempL $(tempL)")
close(fileGro)

#save top file
fileTop=open("test.top","w")
text="""
[ defaults ]
1    2 no    1.00000   1.00000
[ atomtypes ]
O     1.0000    0.0000 A     0.3000000000E+00    0.8314459585E-02
A     1.0000    0.0000 A     0.3000000000E+00    0.8314459585E-02
B     1.0000    0.0000 A     0.3000000000E+00    0.8314459585E-02
C     1.0000    0.0000 A     0.3000000000E+00    0.8314459585E-02
M     1.0000    0.0000 A     0.3000000000E+00    0.8314459585E-02
[ nonbond_params ]
M     A     1     0.3000000000E+00     1.662891917E-02
M     B     1     0.3000000000E+00     0.8314459585E-02

[ moleculetype ]
A    0
[ atoms ]
1 A     1 A  A     1   0.00000   1.00000
[ bonds ]
[ pairs ]
[ angles ]
[ dihedrals ]
[ exclusions ]
[ position_restraints ]
[ constraints ]

[ moleculetype ]
B    0
[ atoms ]
1 B     1 B  B     1   0.00000   1.00000
[ bonds ]
[ pairs ]
[ angles ]
[ dihedrals ]
[ exclusions ]
[ position_restraints ]
[ constraints ]

[ moleculetype ]
C    0
[ atoms ]
1 C     1 C  C     1   0.00000   1.00000
[ bonds ]
[ pairs ]
[ angles ]
[ dihedrals ]
[ exclusions ]
[ position_restraints ]
[ constraints ]

[ moleculetype ]
M    0
[ atoms ]
1 M     1 M  M     1   0.00000   1.00000
[ bonds ]
[ pairs ]
[ angles ]
[ dihedrals ]
[ exclusions ]
[ position_restraints ]
[ constraints ]

[ system ]
generateg with membcreat v5
[ molecules ]
"""

print(fileTop,"$(text)")
if(NAMol>0)
    println(fileTop,"   A  $(NAMol)")
end
if(NBMol>0)
    println(fileTop,"   B  $(NBMol)")
end
println(fileTop,"   M  $(Nmem)")
close(fileTop)
#

#save mdp file
fileMdp=open("test.mdp","w")
text="""
; RUN CONTROL PARAMETERS =
integrator               = steep 
tinit                    = 0.0
dt                       = 0.0005
nsteps                   = 2000
comm-mode                = Linear  ;

; number of steps for center of mass motion removal
nstcomm                  = 10          ; group(s) for center of mass motion removal
;comm-grps               = NTB
;define                  = -DPOSRES
emtol                    = 10.0 
emstep                   = 0.01


; OUTPUT CONTROL OPTIONS =
nstxout                  = 200 ; save trr file
nstvout                  = 200 ;
nstfout                  = 0
nstlog                   = 200
nstenergy                = 200
nstxtcout                = 200 ; xtc file
xtc_precision            = 200

; NEIGHBORSEARCHING PARAMETERS =
nstlist                  = 60      ;как часто обновлять лист соседей
ns_type                  = simple ;grid  ;использовать сетку grid для больших систем, simple проверять все молекулы (дляя маленьких)
pbc                      = xyz ;xyz  - периодические граничные условаия для во все стороны, no  - без периодичности, xy - только для двух сторон  по z - бесконечные
periodic_molecules       = no
rlist                    = 2.0

; OPTIONS FOR ELECTROSTATICS AND VDW =
coulombtype              = Cut-off
rcoulomb                 = 2.0
vdw_type                 = Cut-off
rvdw                     = 2.0 ; nm
fourierspacing           = 0.12 ; nm
DispCorr                 = EnerPres
pme_order                = 4
ewald_rtol               = 1e-05
;ewald_geometry          = 3dc
optimize_fft             = yes
cutoff-scheme            = Verlet


; OPTIONS FOR WEAK COUPLING ALGORITHMS =
tcoupl                   = nose-hoover
"""
text3="""
Pcoupl                   = no 

; GENERATE VELOCITIES FOR STARTUP RUN =
gen_vel                  = no
gen_temp                 = 298.0
gen_seed                 = 473529

; OPTIONS FOR BONDS    
constraints              = none
constraint_algorithm     = lincs 
unconstrained_start      = no
shake_tol                = 0.0001
lincs_order              = 4
lincs_warnangle          = 30
morse                    = no
lincs_iter               = 1

;Non-equilibrium MD
freezegrps		= M
freezedim		= Y Y Y

"""
if((NAMol>0)&&(NBMol>0))
    text2="""tc-grps                  = A B M
tau_t                    = 0.7 0.7  0.7
ref_t                    = 1.5 1.5 1.5
"""
end
if((NAMol>0)&&(NBMol<1))
    text2="""tc-grps                  = A M
tau_t                    = 0.7  0.7
ref_t                    = 1.5 1.5
"""
end
if((NAMol<1)&&(NBMol>0))
    text2="""tc-grps                  = B M
tau_t                    = 0.7  0.7
ref_t                    = 1.5 1.5
"""
end

print(fileMdp,"$(text)")
print(fileMdp,"$(text2)")
print(fileMdp,"$(text3)")
close(fileMdp)
#
fileId=open("start.sh","w")
text="""
rm test.tpr
gmx grompp -f test.mdp -c test.gro -n test.ndx -p test.top -o test.tpr
gmx mdrun -deffnm test 
find -name "test.mdp" -print0 | xargs --null sed -i 's/integrator               = steep/integrator               = md/'
find -name "test.mdp" -print0 | xargs --null sed -i 's/nsteps                   = 2000/nsteps                   = 200000/'
gmx grompp -f test.mdp -c test.gro -n test.ndx -p test.top -o test.tpr
gmx mdrun -deffnm test 
find -name "test.mdp" -print0 | xargs --null sed -i 's/nsteps                   = 200000/nsteps                   = 2000000/'
gmx grompp -f test.mdp -c test.gro -n test.ndx -p test.top -o test.tpr
gmx mdrun -deffnm test 

"""
println(fileId,"$(text)")
close(fileId)

id=0
fileId=open("test.ndx","w")
if(NAMol>0)
    println(fileId," [A]")
    for i=1:NAMol
        global id+=1
        @printf(fileId,"%d  ",id)
        if(mod(id,10)==0)
            @printf(fileId,"\n")
        end
    end
end
if(NBMol>0)
    println(fileId,"\n [B]")
    for i=1:NBMol
        global id+=1
        @printf(fileId,"%d  ",id)
        if(mod(id,10)==0)
            @printf(fileId,"\n")
        end
    end
end
println(fileId,"\n [M]")
for i=1:Nmem
    global id+=1
    @printf(fileId,"%d  ",id)
    if(mod(id,10)==0)
        @printf(fileId,"\n")
    end
end

close(fileId)





