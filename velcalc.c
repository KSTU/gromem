#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define BUFFER 255

int find_maximum(int a[], int n) {
	int c, max, index;
	max = a[0];
	index = 0;
	for (c = 1; c < n; c++) {
		if (a[c] > max) {
			index = c;
			max = a[c];
		}
	}
	return a[index];
}

int get_integer(char a[], int begin, int length) {
	char b[length];
	for (int i = 0; i < length; i++) {
		b[i] = a[begin + i];
	}
	return atoi(b);
}

void get_string(char a[], char out[], int begin, int length) {
	for (int i = 0; i < length; i++) {
		out[i] = a[begin + i];
	}
}

float get_float(char a[], int begin, int length) {
	char b[length];
	for (int i = 0; i < length; i++) {
		b[i] = a[begin + i];
	}
	return atof(b);
}

int main (int argc, char * argv[]){
//varailes
FILE* fileId;
char tempString[BUFFER];
int start;
int aNum;
int chk;
int tempInt;

int* curMol;
int* curAtom;
char** molName;
char** atomName;
float* x;
float* y;
float* z;
float* vx;
float* vy;
float* vz;

float boxX,boxY,boxZ;

char** listMol;
int listNmol;

//molecules array
int* molType;
float* molCX;
float* molCY;
float* molCZ;
float* molCVX;
float* molCVY;
float* molCVZ;

float centerX, centerY, centerZ;
float centerVx, centerVy, centerVz;
int tempInt2, id;

int mNum;

//arrays
int dn0, dn1, dn2, dn3;
float** dens0;
float** dens1;
float** dens2;
float** dens3;
int id0, id1, id2, id3;
float* r0;
float* r1;
float* r2;
float* r3;
float nstep;
float** avvz0;
float** avvz1;
float** avvz2;
float** avvz3;
float** avv0;
float** avv1;
float** avv2;
float** avv3;

	//functions
	int find_maximum(int a[], int n);
	int get_integer(char a[], int begin, int length) ;
	void get_string(char a[], char out[], int begin, int length);
	float get_float(char a[], int begin, int length);

	//read test.temp check new data
	fileId=fopen("test.temp","r");
	for(int i=0;i<16;i++){
		chk=fgets(tempString,BUFFER,fileId);
	}
	chk=fscanf(fileId,"%d",&start);
	if(start==1){
		printf("Start new average\n",start);
	}
	fclose(fileId);
	//read gro file
	if(argc<2){
		printf("usege: velcalc [grofilename]\n");
		return 1;
	}
	fileId=fopen(argv[1],"r");
	if(fileId==NULL){
		printf("Can't open file %s %d\n",argv[1], argc);
		return 1;
	}
	chk=fgets(tempString,BUFFER,fileId);
	chk=fgets(tempString,BUFFER,fileId);
	aNum=atoi(tempString);
//	chk=fscanf(fileId,"%d",&aNum);

	curMol=(int*)malloc(aNum*sizeof(int));
	curAtom=(int*)malloc(aNum*sizeof(int));
	molName=(char**)malloc(aNum*sizeof(char*));
	for(int i=0;i<aNum;i++){
		molName[i]=(char*)malloc(5*sizeof(char));
	}
	atomName=(char**)malloc(aNum*sizeof(char*));
	for(int i=0;i<aNum;i++){
		atomName[i]=(char*)malloc(5*sizeof(char));
	}
	x=(float*)malloc(aNum*sizeof(float));
	y=(float*)malloc(aNum*sizeof(float));
	z=(float*)malloc(aNum*sizeof(float));
	vx=(float*)malloc(aNum*sizeof(float));
	vy=(float*)malloc(aNum*sizeof(float));
	vz=(float*)malloc(aNum*sizeof(float));

	printf("mol number %d\n",aNum);
	for(int i=0;i<aNum;i++){
		fgets(tempString, BUFFER, fileId);    //read string to buffer
/*		printf("%s\n 111",tempString);*/
/*		printf("string %s\n",tempString);*/
		get_string(tempString, molName[i], 5, 5);    //parse string
/*		printf("sub %s\n", molName[i]);*/
		get_string(tempString, atomName[i], 10, 5);
/*		printf("atom %s\n", atomName[i]);*/
		curMol[i]=get_integer(tempString,0,5);
/*		printf("bbb");*/
		curAtom[i]=get_integer(tempString,15,5);
/*		printf("ccc");*/
		x[i] = get_float(tempString, 20, 8);
		y[i] = get_float(tempString, 28, 8);
		z[i] = get_float(tempString, 36, 8);
		vx[i]=get_float(tempString,44,8);
		vy[i]=get_float(tempString,52,8);
		vz[i]=get_float(tempString,60,8);
		printf("mol %d atom %d molname %s atomname %s x %f y %f z %f  \n", curMol[i],curAtom[i], molName[i],atomName[i], vx[i],vy[i],vz[i]);
	}
	fscanf(fileId," %f %f %f",&boxX,&boxY,&boxZ);
	printf("X %f Y %f Z %f \n",boxX,boxY,boxZ);
	fclose(fileId);

//get list of molecule names
	listMol=(char**)malloc(20*sizeof(char*));
	for(int i=0;i<20;i++){
		listMol[i]=(char*)malloc(5*sizeof(char));
	}
	strcpy(listMol[0],molName[0]);
	listNmol=1;
	for(int i=0;i<aNum;i++){
		tempInt=0;
		for(int j=0;j<listNmol;j++){
			if(strcmp(molName[i],listMol[j])==0){
				tempInt=1;
			}
		}
		if(tempInt==0){
			strcpy(listMol[listNmol],molName[i]);
			printf("list %s mol %s\n",listMol[listNmol],molName[i]);
			listNmol++;
		}
	}
	printf("number of substance %d\n",listNmol);
	printf("%s\n",listMol[0]);
	printf("%s\n",listMol[1]);
	printf("%s\n",listMol[2]);
	for(int i=0;i<listNmol;i++){
		printf("%s\n",listMol[i]);
	}
//get center of molecules
	mNum=curMol[aNum-1];
	molType=(int*)malloc(mNum*sizeof(int));
	molCX=(float*)malloc(mNum*sizeof(float));
	molCY=(float*)malloc(mNum*sizeof(float));
	molCZ=(float*)malloc(mNum*sizeof(float));
	molCVX=(float*)malloc(mNum*sizeof(float));
	molCVY=(float*)malloc(mNum*sizeof(float));
	molCVZ=(float*)malloc(mNum*sizeof(float));
	
	//current average
	tempInt=curMol[0];
	tempInt2=1;
	centerX=x[0];
	centerY=y[0];
	centerZ=z[0];
	centerVx=vx[0];
	centerVy=vy[0];
	centerVz=vz[0];
	id=0;
	for(int i=1;i<aNum;i++){
		if(curMol[i]!=tempInt){	//если молекула другая то
			tempInt++;
			molCX[id]=centerX/tempInt2;
			molCY[id]=centerY/tempInt2;
			molCZ[id]=centerZ/tempInt2;
			molCVX[id]=centerVx/tempInt2;
			molCVY[id]=centerVy/tempInt2;
			molCVZ[id]=centerVz/tempInt2;
			for(int j=0;j<listNmol;j++){
				if(strcmp(molName[i],listMol[j])==0){
					molType[id]=j;
				}
			}
			tempInt2=0;
			id++;
			centerX=0;
			centerY=0;
			centerZ=0;
			centerVx=0;
			centerVy=0;
			centerVz=0;
		}
		centerX+=x[i];
		centerY+=y[i];
		centerZ+=z[i];
		centerVx+=vx[i];
		centerVy+=vy[i];
		centerVz+=vz[i];
		tempInt2++;
	}

	

//усреднение
	dn0=50; 
	dn1=100;
	dn2=200;
	dn3=300;

	//debug output to file
	fileId=fopen("mol.out","w");
		for(int i=0;i<curMol[aNum-1];i++){
			fprintf(fileId,"%d\t%f\t%f\t%f\n",molType[i],molCZ[i]/boxZ*dn0,molCVY[i],molCVZ[i]);
		}
	fclose(fileId);

	printf("test\n");
	dens0=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		dens0[i]=(float*)calloc(dn0,sizeof(float));
	}
	
	dens1=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		dens1[i]=(float*)calloc(dn1,sizeof(float));
	}
	dens2=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		dens2[i]=(float*)calloc(dn2,sizeof(float));
	}
	dens3=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		dens3[i]=(float*)calloc(dn3,sizeof(float));
	}
	r0=(float*)malloc(dn0*sizeof(float));
	r1=(float*)malloc(dn1*sizeof(float));
	r2=(float*)malloc(dn2*sizeof(float));
	r3=(float*)malloc(dn3*sizeof(float));

	//z velosity
	avvz0=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avvz0[i]=(float*)calloc(dn0,sizeof(float));
	}
	avvz1=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avvz1[i]=(float*)calloc(dn1,sizeof(float));
	}
	avvz2=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avvz2[i]=(float*)calloc(dn2,sizeof(float));
	}
	avvz3=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avvz3[i]=(float*)calloc(dn3,sizeof(float));
	}
	//v^2 averages
	avv0=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avv0[i]=(float*)calloc(dn0,sizeof(float));
	}
	avv1=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avv1[i]=(float*)calloc(dn1,sizeof(float));
	}
	avv2=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avv2[i]=(float*)calloc(dn2,sizeof(float));
	}
	avv3=(float**)malloc(listNmol*sizeof(float*));
	for(int i=0;i<listNmol;i++){
		avv3[i]=(float*)calloc(dn3,sizeof(float));
	}

//
	
	
	nstep=0.0;
	
	for(int i=0;i<dn0;i++){
		r0[i]=(i+0.5)/dn0*boxZ;
	}
	for(int i=0;i<dn1;i++){
		r1[i]=(i+0.5)/dn1*boxZ;
	}
	for(int i=0;i<dn2;i++){
		r2[i]=(i+0.5)/dn2*boxZ;
	}
	for(int i=0;i<dn3;i++){
		r3[i]=(i+0.5)/dn3*boxZ;
	}
	
	if(start!=1){
		//read from file
		fileId=fopen("averages.temp","r");
		fscanf(fileId,"%f\n",&nstep);
		for(int i=0;i<dn0;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&dens0[j][i]);
			}
		}
		for(int i=0;i<dn1;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&dens1[j][i]);
			}
		}
		for(int i=0;i<dn2;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&dens2[j][i]);
			}
		}
		for(int i=0;i<dn3;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&dens3[j][i]);
			}
		}
		for(int i=0;i<dn0;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avvz0[j][i]);
			}
		}
		for(int i=0;i<dn1;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avvz1[j][i]);
			}
		}
		for(int i=0;i<dn2;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avvz2[j][i]);
			}
		}
		for(int i=0;i<dn3;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avvz3[j][i]);
			}
		}
		///
		for(int i=0;i<dn0;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avv0[j][i]);
			}
		}
		for(int i=0;i<dn1;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avv1[j][i]);
			}
		}
		for(int i=0;i<dn2;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avv2[j][i]);
			}
		}
		for(int i=0;i<dn3;i++){
			for(int j=0;j<listNmol;j++){
				fscanf(fileId,"%f",&avv3[j][i]);
			}
		}
		fclose(fileId);
	}
	printf("test2\n");
	for(int i=0;i<mNum;i++){
		id0=floor(molCZ[i]/boxZ*dn0);
		printf("%d %d %f\n",i,id0,molCZ[i]/boxZ*dn0);
		id1=floor(molCZ[i]/boxZ*dn1);
		id2=floor(molCZ[i]/boxZ*dn2);
		id3=floor(molCZ[i]/boxZ*dn3);

		dens0[molType[i]][id0]+=1.0;
		dens1[molType[i]][id1]+=1.0;
		dens2[molType[i]][id2]+=1.0;
		dens3[molType[i]][id3]+=1.0;

		avvz0[molType[i]][id0]+=molCVZ[i];
		avvz1[molType[i]][id1]+=molCVZ[i];
		avvz2[molType[i]][id2]+=molCVZ[i];
		avvz3[molType[i]][id3]+=molCVZ[i];

		avv0[molType[i]][id0]+=molCVX[i]*molCVX[i]+molCVY[i]*molCVY[i]+molCVZ[i]*molCVZ[i];
		avv1[molType[i]][id1]+=molCVX[i]*molCVX[i]+molCVY[i]*molCVY[i]+molCVZ[i]*molCVZ[i];
		avv2[molType[i]][id2]+=molCVX[i]*molCVX[i]+molCVY[i]*molCVY[i]+molCVZ[i]*molCVZ[i];
		avv3[molType[i]][id3]+=molCVX[i]*molCVX[i]+molCVY[i]*molCVY[i]+molCVZ[i]*molCVZ[i];
	}
	nstep+=1.0;
	printf("test3\n");
	
//plot averages
	fileId=fopen("v-dens0.out","w");
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn0;i++){
		fprintf(fileId,"%f",r0[i]);
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"\t%f",dens0[j][i]/nstep/(boxX*boxY*boxZ/dn0));
		}
		fprintf(fileId,"\n");
	}
	fclose(fileId);
	//1
	fileId=fopen("v-dens1.out","w");
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn1;i++){
		fprintf(fileId,"%f",r1[i]);
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"\t%f",dens1[j][i]/nstep/(boxX*boxY*boxZ/dn1));
		}
		fprintf(fileId,"\n");
	}
	fclose(fileId);
	//2
	fileId=fopen("v-dens2.out","w");
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn2;i++){
		fprintf(fileId,"%f",r2[i]);
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"\t%f",dens2[j][i]/nstep/(boxX*boxY*boxZ/dn2));
		}
		fprintf(fileId,"\n");
	}
	fclose(fileId);
	//3
	fileId=fopen("v-dens3.out","w");
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn3;i++){
		fprintf(fileId,"%f",r3[i]);
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"\t%f",dens3[j][i]/nstep/(boxX*boxY*boxZ/dn3));
		}
		fprintf(fileId,"\n");
	}
	fclose(fileId);
//z-velosity averages
	fileId=fopen("v-zvel0.out","w");	//0
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn0;i++){
		fprintf(fileId,"%f",r0[i]);
		for(int j=0;j<listNmol;j++){
			if(dens0[j][i]>0.0){
				fprintf(fileId,"\t%f",avvz0[j][i]/dens0[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-zvel1.out","w");	//1
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn1;i++){
		fprintf(fileId,"%f",r1[i]);
		for(int j=0;j<listNmol;j++){
			if(dens1[j][i]>0.0){
				fprintf(fileId,"\t%f",avvz1[j][i]/dens1[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-zvel2.out","w");	//2
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn2;i++){
		fprintf(fileId,"%f",r2[i]);
		for(int j=0;j<listNmol;j++){
			if(dens2[j][i]>0.0){
				fprintf(fileId,"\t%f",avvz2[j][i]/dens2[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-zvel3.out","w");	//3
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn3;i++){
		fprintf(fileId,"%f",r3[i]);
		for(int j=0;j<listNmol;j++){
			if(dens3[j][i]>0.0){
				fprintf(fileId,"\t%f",avvz3[j][i]/dens3[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}

	fclose(fileId);

//velosity2 averages
	fileId=fopen("v-vel0.out","w");	//0
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn0;i++){
		fprintf(fileId,"%f",r0[i]);
		for(int j=0;j<listNmol;j++){
			if(dens0[j][i]>0.0){
				fprintf(fileId,"\t%f",avv0[j][i]/dens0[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-vel1.out","w");	//1
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn1;i++){
		fprintf(fileId,"%f",r1[i]);
		for(int j=0;j<listNmol;j++){
			if(dens1[j][i]>0.0){
				fprintf(fileId,"\t%f",avv1[j][i]/dens1[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-vel2.out","w");	//2
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn2;i++){
		fprintf(fileId,"%f",r2[i]);
		for(int j=0;j<listNmol;j++){
			if(dens2[j][i]>0.0){
				fprintf(fileId,"\t%f",avv2[j][i]/dens2[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}
	fileId=fopen("v-vel3.out","w");	//3
	fprintf(fileId,"r");
	for(int i=0;i<listNmol;i++){
		fprintf(fileId,"\t%s",listMol[i]);
	}
	fprintf(fileId,"\n");
	for(int i=0;i<dn3;i++){
		fprintf(fileId,"%f",r3[i]);
		for(int j=0;j<listNmol;j++){
			if(dens3[j][i]>0.0){
				fprintf(fileId,"\t%f",avv3[j][i]/dens3[j][i]);
			} else{
				fprintf(fileId,"\t%f",0.0);
			}
		}
		fprintf(fileId,"\n");
	}

	fclose(fileId);
	
//save averages
	fileId=fopen("averages.temp","w");
	fprintf(fileId,"%f\n",nstep);
	for(int i=0;i<dn0;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",dens0[j][i]);
		}
	}
	for(int i=0;i<dn1;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",dens1[j][i]);
		}
	}
	for(int i=0;i<dn2;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",dens2[j][i]);
		}
	}
	for(int i=0;i<dn3;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",dens3[j][i]);
		}
	}
	//z vel
	for(int i=0;i<dn0;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avvz0[j][i]);
		}
	}
	for(int i=0;i<dn1;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avvz1[j][i]);
		}
	}
	for(int i=0;i<dn2;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avvz2[j][i]);
		}
	}
	for(int i=0;i<dn3;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avvz3[j][i]);
		}
	}
	//v^2
	for(int i=0;i<dn0;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avv0[j][i]);
		}
	}
	for(int i=0;i<dn1;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avv1[j][i]);
		}
	}
	for(int i=0;i<dn2;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avv2[j][i]);
		}
	}
	for(int i=0;i<dn3;i++){
		for(int j=0;j<listNmol;j++){
			fprintf(fileId,"%f\n",avv3[j][i]);
		}
	}
	
	fclose(fileId);
	
	return 0;
}
