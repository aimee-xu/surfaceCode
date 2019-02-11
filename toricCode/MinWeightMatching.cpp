//
//  example of accessing blossom code for Minimum Weight Perfect Matching
//
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "PerfectMatching.h"
#include "allFunctions.h"



//--------------------------------------------------------------
//---------------------- START OF main( )  ----------------------
//--------------------------------------------------------------
int main ( int argc, char *argv[] )
{
    
    time_t rSeed;
    rSeed = time (NULL);
    unsigned long labelNumber=time (NULL)*10000000;
    srand(labelNumber);
    srand( rand()+ rSeed);
    
    double errorRateGate=0.008;
    double errorRateMeasure=0.008;
    int index;
    
    
//==========================================================
    //generate look-up table
    
    int patternX[10]={0};
    int patternZ[10]={0};
    
    long int *storeX; long int numofstore = 1L << 10;
    storeX=(long int *)calloc(numofstore,sizeof(long int));
    if (storeX == NULL)
    { printf("space failed to be allocated！2\n"); exit(1);}
    
    long int *storeZ;
    storeZ=(long int *)calloc(numofstore,sizeof(long int));
    if (storeZ == NULL)
    { printf("space failed to be allocated！2\n"); exit(1);}
    
    int countCycle=0; long int lookup=1000000;
    while (countCycle<lookup){
        
        measureNoise(4,patternZ,errorRateMeasure);
        cnot_r (0,4, patternZ, errorRateGate);
        cnot_r (1,4, patternZ, errorRateGate);
        cnot_r (2,4, patternZ, errorRateGate);
        cnot_r (3,4, patternZ, errorRateGate);
        measureNoise(4,patternZ,errorRateMeasure);
        
        int temp=0;
        
        for(index=8;index>=0;index--)
        {
            int coef= 1L << index;
            temp+=coef*patternZ[index];
            
        }
        storeZ[temp]++;
        
        measureAbsolute_r(4,patternZ);
        measureAbsolute_r(3,patternZ);
        measureAbsolute_r(2,patternZ);
        measureAbsolute_r(1,patternZ);
        measureAbsolute_r(0,patternZ);
        
        countCycle++;
        
    }
    
    countCycle=0;
    while (countCycle<lookup){
        
        measureNoise(4,patternX, errorRateMeasure);
        hadamard_r (4, patternX, errorRateGate);
        cnot_r (4,0, patternX, errorRateGate);
        cnot_r (4,1, patternX, errorRateGate);
        cnot_r (4,2, patternX, errorRateGate);
        cnot_r (4,3, patternX, errorRateGate);
        hadamard_r (4, patternX, errorRateGate);
        measureNoise(4,patternX, errorRateMeasure);
        
        int temp=0;
        
        for(index=8;index>=0;index--)
        {
            int coef= 1L << index;
            temp+=coef*patternX[index];
            
        }
        storeX[temp]++;
        
        measureAbsolute_r(4,patternX);
        measureAbsolute_r(3,patternX);
        measureAbsolute_r(2,patternX);
        measureAbsolute_r(1,patternX);
        measureAbsolute_r(0,patternX);
        
        countCycle++;
        
    }
   
    
    //===========================================================
    //main simulation
    long int modelx=0; long int modelz=0;
    int iterations=13; int numCycle=10000;
    int numQ=64;
    int n= sqrt(numQ);
    int i; int j; int p; int q; int t; int cycle;
    
    //space allocation
    
    int *globalStore;
    globalStore=(int *)calloc(numQ*2,sizeof(int));
    if (globalStore == NULL){
        printf("space failed to be allocated！2\n"); exit(1);}
    
    int *copyX;
    copyX=(int *)calloc(numQ*2,sizeof(int));
    if (copyX == NULL){
        printf("space failed to be allocated！2\n"); exit(1);}
    
    int *copyZ;
    copyZ=(int *)calloc(numQ*2,sizeof(int));
    if (copyZ == NULL){
        printf("space failed to be allocated！2\n"); exit(1);}
    
    
    int weightTime=1;
    int weightSpace=1;
    
    
    int rightX=0; int rightZ=0;
    
//-----------------------------------------------------

    for (cycle=0;cycle<numCycle;cycle++){
        
        
        int labelTableX[2000]={0};
        int labelTableZ[2000]={0};
        
        int tableX[1000][1000]={0};
        int tableZ[1000][1000]={0};
        
        int labelX=0;
        int labelZ=0;
        
        for (p=0;p<2*numQ;p++){
            globalStore[p]=0;
            copyX[p]=0;
            copyZ[p]=0;
        }
            
        int count=1;
        
    while (count<iterations){

        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                t=2*n*p+2*q;
                check (t, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, storeZ,lookup,1);
            }
            p++;
        }
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                t=n+1+2*n*p+2*q;
                check (t, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, storeX,lookup,1);
            }
            p++;
        }

        
        count++;
    }

    //the last cycle
        
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2;q++){
                t=2*n*p+2*q;
                check (t, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, storeZ, lookup,0);
            }
            p++;
        }
        
    q=0; p=0;
    while (p<(n/2)){
        for (q=0;q<n/2;q++){
            t=n+1+2*n*p+2*q;
            check (t, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, storeX, lookup,0);
        }
        p++;
    }
        
   
    
    // table x
    creatTable (labelTableX, tableX, weightTime, weightSpace,numQ);
    // table z
    creatTable (labelTableZ, tableZ, weightTime, weightSpace,numQ);
        
//---------------------------------------------------------------------
        
    int numNodesX=labelX;
    int numNodesZ=labelZ;
        modelx+=labelX;
        modelz+=labelZ;
    int numEdgesX=0;
    int numEdgesZ=0;
        
    for (p=0;p<labelX;p++){
        for (q=p+1;q<labelX;q++){
            if (tableX[p][q]!=0)
                numEdgesX++;
        }
    }
        
    for (p=0;p<labelZ;p++){
        for (q=p+1;q<labelZ;q++){
            if (tableZ[p][q]!=0)
                numEdgesZ++;
        }
    }
    
    int workEdgeX=numEdgesX*2;
    int workEdgeZ=numEdgesZ*2;

//        for (p=0;p<n;p++){
//            for (q=0;q<n;q++)
//                printf("%d ",globalStore[2*(p*n+q)+1]);
//            printf("\n");
//        }
//        printf("\n");
    
    //the lines below initialise the system for perfect matching
	struct PerfectMatching::Options options;
	PerfectMatching *pmX;
	pmX= new PerfectMatching(numNodesX, workEdgeX);
	
    //the lines below switch off the information that is otherwise printed out by the perfect matching system
    options.verbose = false;
    pmX->options = options;
    
    //now we just load the graph information into the system
    for (i=0; i<numNodesX-1; i++)
    	for (j=i+1; j<numNodesX; j++)
    		pmX->AddEdge(i, j, tableX[i][j]);
    		
	pmX->Solve(); //this instruction tell the perfect matching system to find the solution
        
	for (i=0; i<numNodesX; i++){
		int j = pmX->GetMatch(i); //this line just gets the partner for node i from the system
        if (j>i)
            correctionX(i,j,n,labelTableX,globalStore);
    }
//        sortError(globalStore,numQ);
//        for (p=0;p<n;p++){
//            for (q=0;q<n;q++)
//                printf("%d ",globalStore[2*(p*n+q)+1]);
//            printf("\n");
//        }
//        printf("\n");
    
    for (p=0;p<2*numQ;p++)
        globalStore[p]=globalStore[p]%2;
        
    int result=sortLogicalZ(globalStore,numQ);
    if (result==0)
        rightX++;
        
//----------------------------------------------------------
        
//        struct PerfectMatching::Options options;
        PerfectMatching *pm;
        pm= new PerfectMatching(numNodesZ, workEdgeZ);
        
        //the lines below switch off the information that is otherwise printed out by the perfect matching system
        options.verbose = false;
        pm->options = options;
        
        //now we just load the graph information into the system
        for (i=0; i<numNodesZ-1; i++)
            for (j=i+1; j<numNodesZ; j++)
                pm->AddEdge(i, j, tableZ[i][j]);
        
        pm->Solve(); //this instruction tell the perfect matching system to find the solution
        
        for (i=0; i<numNodesZ; i++){
            int j = pm->GetMatch(i); //this line just gets the partner for node i from the system
            if (j>i)
                correctionZ(i,j,n,labelTableZ,globalStore);
        }
        
        
        for (p=0;p<2*numQ;p++)
            globalStore[p]=globalStore[p]%2;
        
        result=sortLogicalX(globalStore,numQ);
        if (result==0)
            rightZ++;

    }
    
//    printf("model=%ld\n",model);
    double ratio1=(double)rightX/(double)numCycle;
    double ratio2=(double)rightZ/(double)numCycle;
    printf("%f\n",ratio1);
    printf("%f\n",ratio2);
    printf("modelx=%ld\n",modelx);
    printf("modelz=%ld\n",modelz);

  
    
    free(globalStore);
    free(copyX);
    free(copyZ);
    free(storeX);
    free(storeZ);
    
    return 0;
}
