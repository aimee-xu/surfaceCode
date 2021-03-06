#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "PerfectMatching.h"
#include "allFunctionsPlanerS.h"

//--------------------------------------------------------------
//---------------------- START OF main( )  ----------------------
//--------------------------------------------------------------
int main ( int argc, char *argv[] )
{
    
    time_t rSeed;
    rSeed = time (NULL);
    unsigned long labelNumber=0;
    
    int iterations=15;
    int numQ=64;
    int width=2; int numCycle=1000;
    double gateError=0;
    int index; int jndex; int p; int q; int ts;

    // begin paste -------------------------------
    if (argc<3){ //this message will print if there is not at least one command line parameter
        printf("\nINFORMATION ON COMMAND LINE PARAMETERS:\nA usage example is './thisProg -r=5 -i=200 -em=0.1 -eg=0.1 -ep=0.1 -t=200'\n");
        printf("  -i is the numcycle\n");
        printf("  -c is the number of iterations\n");
        printf("  -a is errorrate ancillar\n");
        printf("  -d is error rate data qubit\n");
        printf("  -u is number of qubit\n");
        exit(1);
    }
    
    int i;
    for (i=1; i<argc;  i++){
        if (argv[i][0]=='-'){
            int j=1;
            while (argv[i][j]=='-') j++;
            int failedToParse=1;
            char thisFlag=argv[i][j]; j++;
            if (argv[i][j]=='='){
                while (j>=0){ 		//this loop erases the '-b=' in a string like '-b=1.3' so that only '   1.3' remains, so it can be converted to a number cleanly
                    argv[i][j]=' ';
                    j--;
                }
                if (thisFlag=='c'){
                    iterations=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='i'){
                    numCycle=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='u'){
                    numQ=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='d'){
                    width=atoi(argv[i]);
                    failedToParse=0;
                }
                if (thisFlag=='g'){
                    gateError=strtof(argv[i], NULL);
                    failedToParse=0;
                }
                if (thisFlag=='r'){
                    labelNumber=atoi(argv[i])*1000;
                    failedToParse=0;
                }
            }
            if (failedToParse==1){ printf("\nCannot interpret command line argument %s. (Note use '=' as in '-s=100'). ABORTING.\n",argv[i]); exit(1); }
            
        }else{ // an argument not starting with - detected
            printf("A command line parameter %s seen; cannot interpret this. ABORTING.\n",argv[i]); exit(1);
        }
    }
    printf(" The settings are -i=%ld -c=%d -r=%d -u=%d\"; \n", numCycle,iterations,labelNumber,numQ);
    srand(labelNumber);
    srand( rand()+ rSeed);
    
//==========================================================
    //load data
    double CnotPro[400000]={0};

    char nameWidth[21];
    char tempi[14]=".long12.3.txt";
    sprintf(nameWidth, "%d", width);
    strcat(nameWidth,tempi);

    
    FILE *fp;
    if ((fp=fopen(nameWidth,"r"))==NULL){ printf("cannot open file\n"); exit(1);}
    for (index=0;index<400000;index++)
        fscanf(fp,"%lf",&CnotPro[index]);
    fclose(fp);
    
//==========================================================
    //main simulation
    int n= sqrt(numQ);
    int cycle;
    
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
    
    int right=0;
    //-----------------------------------------------------
    for (cycle=0;cycle<numCycle;cycle++){
        
        int labelTableX[6000]={0};
        int labelTableZ[6000]={0};
        int labelX=0;
        int labelZ=0;
        
        for (p=0;p<2*numQ;p++){
            globalStore[p]=0;
            copyX[p]=0;
            copyZ[p]=0;
        }
        
        int count=1;
        
        int *store;
        store=(int *)calloc(numQ,sizeof(int));
        for (index=0; index<numQ; index++){
            int ran =(rand() % 50000);
            store[index]=ran;
        }
        
        while (count<iterations){
            
            q=0; p=0;
            while (p<(n/2)-1){
                for (q=0;q<n/2;q++){
                    ts=n+1+2*n*p+2*q;
                    if (q==0)
                        stabiliser (ts, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, store,CnotPro,1,0,gateError);//1st0:x; 2nd1:add error 3st0:direction
                    else if (q==n/2-1)
                        stabiliser (ts, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, store,CnotPro,1,1,gateError);
                    else
                        stabiliser (ts, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, store,CnotPro,1,3,gateError);
                }
                p++;
            }
            
            q=0; p=0;
            while (p<(n/2)){
                for (q=0;q<n/2-1;q++){
                    ts=2*n*p+2*q;
                    if (p==0)
                        stabiliser (ts, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, store,CnotPro,1,0,gateError);
                    else if (p==n/2-1)
                        stabiliser (ts, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, store,CnotPro,1,1,gateError);
                    else
                        stabiliser (ts, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, store,CnotPro,1,3,gateError);
                }
                p++;
            }
            
            count++;
        }
        
        //the last cycle
        q=0; p=0;
        while (p<(n/2)-1){
            for (q=0;q<n/2;q++){
                ts=n+1+2*n*p+2*q;
                stabiliser (ts, 0, copyX, globalStore, count, numQ, &labelX, labelTableX, store,CnotPro,0,3,gateError);
            }
            p++;
        }
        
        q=0; p=0;
        while (p<(n/2)){
            for (q=0;q<n/2-1;q++){
                ts=2*n*p+2*q;
                stabiliser (ts, 1, copyZ, globalStore, count, numQ, &labelZ, labelTableZ, store,CnotPro,0,3,gateError);
            }
            p++;
        }
        
        int *tableX;
        tableX=(int *)calloc(labelX*labelX*3,sizeof(int));
        if (tableX == NULL){
            printf("space failed to be allocated！2\n"); exit(1);}
        
        int *tableZ;
        tableZ=(int *)calloc(labelZ*labelZ*3,sizeof(int));
        if (tableZ == NULL){
            printf("space failed to be allocated！2\n"); exit(1);}
        
        // table x
        int numEdgesX=creatTable (labelTableX,labelX,tableX,weightTime,weightSpace,numQ,0);
        // table z
        int numEdgesZ=creatTable (labelTableZ,labelZ,tableZ,weightTime,weightSpace,numQ,1);
        
        //---------------------------------------------------------------------
        
        int numNodesX=2*labelX;
        int numNodesZ=2*labelZ;
        
        if (numEdgesX!=labelX*labelX){printf("numEdgeX wrong\n");exit(1);}
        if (numEdgesZ!=labelZ*labelZ){printf("numEdgeZ wrong\n");exit(1);}
        
        int workEdgeX=numEdgesX*2;
        int workEdgeZ=numEdgesZ*2;
        
        //the lines below initialise the system for perfect matching
        struct PerfectMatching::Options options;
        PerfectMatching *pmX;
        pmX= new PerfectMatching(numNodesX, workEdgeX);
        
        //the lines below switch off the information that is otherwise printed out by the perfect matching system
        options.verbose = false;
        pmX->options = options;
        
        //now we just load the graph information into the system
        for (index=0; index<numEdgesX; index++)
            pmX->AddEdge(tableX[3*index], tableX[3*index+1], tableX[3*index+2]);
        
        pmX->Solve(); //this instruction tell the perfect matching system to find the solution
        for (index=0; index<numNodesX; index++){
            int jndex = pmX->GetMatch(index); //this line just gets the partner for node i from the system
            if (jndex>index)
                errorCorrect(n,jndex,index,labelX,labelTableX,0,globalStore);
        }
        
        for (p=0;p<2*numQ;p++)
            globalStore[p]=globalStore[p]%2;
        
        int resultx=sortLogicalZ(globalStore,numQ);
        if (resultx==0)
            rightX++;
        
        //--------------------------------------------------------
        PerfectMatching *pm;
        pm= new PerfectMatching(numNodesZ, workEdgeZ);
        
        //the lines below switch off the information that is otherwise printed out by the perfect matching system
        options.verbose = false;
        pm->options = options;
        
        //now we just load the graph information into the system
        for (index=0; index<numEdgesZ; index++)
            pm->AddEdge(tableZ[3*index], tableZ[3*index+1], tableZ[3*index+2]);
        
        pm->Solve(); //this instruction tell the perfect matching system to find the solution
        for (index=0; index<numNodesZ; index++){
            int jndex = pm->GetMatch(index); //this line just gets the partner for node i from the system
            if (jndex>index)
                errorCorrect(n,jndex,index,labelZ,labelTableZ,1,globalStore);
        }
    
        for (p=0;p<2*numQ;p++)
            globalStore[p]=globalStore[p]%2;
        
        int resultz=sortLogicalX(globalStore,numQ);
//        if (result==0)
//            rightZ++;
        
        if (resultz+resultx==0)
            right++;
        
        delete pm;
        delete pmX;
        
        free(tableZ);
        free(tableX);
        free(store);
        
    }
    
    double ratio=(double)(right)/(double)numCycle;
    printf("%f\n",ratio);
    
    
    free(globalStore);
    free(copyX);
    free(copyZ);

    return 0;
}
