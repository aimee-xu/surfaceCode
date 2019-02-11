#ifndef _ALLFUNCTIONS_H_
#define _ALLFUNCTIONS_H_

# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <complex.h>

const long double oneOverRoot2= 0.7071067811865475244008443621048490392848;
const long double Pi = 3.14159265358979323846264338327950288419716939937510;
const long double oneoverthree=0.3333333333333333333333333333333333333333;

void sortCnotError(int q1, int q2, int *Q)
{
    if (*(Q+2*q1)==1)
        *(Q+2*q2)+=1;
    if (*(Q+2*q2+1)==1)
        *(Q+2*q1+1)+=1;
    *(Q+2*q2)=*(Q+2*q2)%2; *(Q+2*q1+1)=*(Q+2*q1+1)%2;
}

void sortHadamardError(int q, int *Q)
{
    if (*(Q+2*q)==1 && *(Q+2*q+1)==0)
    { *(Q+2*q+1)=1; *(Q+2*q)=0;}
    else if (*(Q+2*q)==0 && *(Q+2*q+1)==1)
    { *(Q+2*q+1)=0; *(Q+2*q)=1;}
}

void measureNoise(int q, int *Q, double gateError)
{
    double ran1; double ran2;
    ran1 =((double)rand()/(double)RAND_MAX);
    ran2 =((double)rand()/(double)RAND_MAX);
    
    if (ran1<gateError){
        if (ran2<=oneoverthree)
        { *(Q+2*q)+=1;}
        else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
        {*(Q+2*q)+=1;*(Q+2*q+1)+=1;}
        else if (ran2>2*oneoverthree)
        {*(Q+2*q+1)+=1;}
        
        *(Q+2*q+1)=*(Q+2*q+1)%2; *(Q+2*q)=*(Q+2*q)%2;
    }
}



void hadamard_r(int q, int *Q, double gateError)
{
    
    sortHadamardError(q,Q);
    double ran1; double ran2;
    ran1 =((double)rand()/(double)RAND_MAX);
    ran2 =((double)rand()/(double)RAND_MAX);
    
    if (ran1<gateError){
        if (ran2<=oneoverthree)
        { *(Q+2*q)+=1;}
        else if (ran2>oneoverthree && ran2<= 2*oneoverthree)
        {*(Q+2*q)+=1;*(Q+2*q+1)+=1;}
        else if (ran2>2*oneoverthree)
        {*(Q+2*q+1)+=1;}
        
        *(Q+2*q+1)=*(Q+2*q+1)%2; *(Q+2*q)=*(Q+2*q)%2;
    }
}

void cnot_r(int q1, int q2, int *Q, double gateError)
{
    
    sortCnotError(q1,q2,Q);
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    if (ran < gateError)
    {
        int myRnd=1 + rand()%15;
        int firstError=myRnd/4;
        int secndError=myRnd%4;
        
        if (firstError==1) {*(Q+2*q1)+=1;}
        if (firstError==2) {*(Q+2*q1)+=1; *(Q+2*q1+1)+=1;}
        if (firstError==3) {*(Q+2*q1+1)+=1;}
        if (secndError==1) {*(Q+2*q2)+=1;}
        if (secndError==2) {*(Q+2*q2)+=1; *(Q+2*q2+1)+=1;}
        if (secndError==3) {*(Q+2*q2+1)+=1;}
        
        *(Q+2*q1)=*(Q+2*q1)%2; *(Q+2*q1+1)=*(Q+2*q1+1)%2; *(Q+2*q2)=*(Q+2*q2)%2; *(Q+2*q2+1)=*(Q+2*q2+1)%2;
        
    }
    
}
void sortNoError (int numQubits, int q, int *Q, double __complex__ *amp)
{
    if (*(Q+2*q)==1 && *(Q+2*q+2)==1){
        *(Q+2*q)=0; *(Q+2*q+2)=0;}
    if (*(Q+2*q+4)==1 && *(Q+2*q+6)==1){
        *(Q+2*q+4)=0; *(Q+2*q+6)=0;}
    if (*(Q+2*q+1)==1 && *(Q+2*q+3)==1 && *(Q+2*q+5)==1 && *(Q+2*q+7)==1){
        *(Q+2*q+1)=0; *(Q+2*q+3)=0;  *(Q+2*q+5)=0;  *(Q+2*q+7)=0;}
    
}

void measureAbsolute_r (int q, int *Q)
{
    *(Q+2*q)=0; *(Q+2*q+1)=0;
}


void measurePerfectX (int *globalStore, int q4, int numQ)
{
    int q1; int q2; int q3; int qa = q4-1;
    int n = sqrt(numQ);
    q1 = q4-n-1;
    q2 = q4+n-1;
    if (q2>=numQ) q2=q2-numQ;
    if (q4%(2*n)==(n+1))
        q3=q4+n-2;
    else
        q3=q4-2;
    
    *(globalStore+2*qa)=0;
    *(globalStore+2*qa+1)=0;
    
    sortCnotError (qa, q1, globalStore);
    sortCnotError (qa, q2, globalStore);
    sortCnotError (qa, q3, globalStore);
    sortCnotError (qa, q4, globalStore);
    
}

void measurePerfectZ (int *globalStore, int q3, int numQ)
{
    int q1; int q2; int q4; int qa = q3+1;
    int n = sqrt(numQ);
    q1 = q3-n+1;
    if (q1<0) q1+=numQ;
    q2 = q3+n+1;
    if (q3%(2*n)==(n-2))
        q4=q3-n+2;
    else
        q4=q3+2;
    
    *(globalStore+2*qa)=0;
    *(globalStore+2*qa+1)=0;
    
    sortCnotError (q1, qa, globalStore);
    sortCnotError (q2, qa, globalStore);
    sortCnotError (q3, qa, globalStore);
    sortCnotError (q4, qa, globalStore);
}


void sortError(int *globalStore, int numQ)
{
    int i; int j; int n =sqrt(numQ);
    for (i=0;i<n/2;i++){
        for (j=0;j<n/2;j++){
            *(globalStore+(1+2*i*n+2*j)*2)=0;
            *(globalStore+(1+2*i*n+2*j)*2+1)=0;
            *(globalStore+(n+2*i*n+2*j)*2)=0;
            *(globalStore+(n+2*i*n+2*j)*2+1)=0;
        }
    }
    
}

void updateMeasure(int *temp, int q, int type, int numQ, int *globalStore)
{
    int q1; int q2; int q3; int q4;int qa; int n = sqrt(numQ);
    
    if (type==0){
        
        q4 = q;
        qa = q-1;
        q1 = q-n-1;
        q2 = q+n-1;
        if (q2>=numQ) q2=q2-numQ;
        if (q%(2*n)==(n+1))
            q3=q+n-2;
        else
            q3=q-2;
    }
    
    if (type==1){

        q3 = q;
        qa = q+1;
        q1 = q-n+1;
        if (q1<0) q1+=numQ;
        q2 = q3+n+1;
        if (q3%(2*n)==(n-2))
            q4=q3-n+2;
        else
            q4=q3+2;
    }
    
    double ran;
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
        *(globalStore+2*q1)+=temp[0]; *(globalStore+2*q1)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q1+1)+=temp[1]; *(globalStore+2*q1+1)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q2)+=temp[2]; *(globalStore+2*q2)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q2+1)+=temp[3]; *(globalStore+2*q2+1)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q3)+=temp[4]; *(globalStore+2*q3)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q3+1)+=temp[5]; *(globalStore+2*q3+1)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q4)+=temp[6]; *(globalStore+2*q4)%=2;
    }
    ran=((double)rand()/(double)RAND_MAX);
    if (ran<2*oneoverthree){
    *(globalStore+2*q4+1)+=temp[7]; *(globalStore+2*q4+1)%=2;
    }
    
    if (type==1){
        ran=((double)rand()/(double)RAND_MAX);
        if (ran<2*oneoverthree){
            *(globalStore+2*qa)+=temp[8]; *(globalStore+2*qa)%=2;
        }
    }
    if (type==0){
        ran=((double)rand()/(double)RAND_MAX);
        if (ran<2*oneoverthree){
            *(globalStore+2*qa+1)+=temp[8]; *(globalStore+2*qa+1)%=2;
        }
    }
}

void check (int q, int flag, int *copy, int *globalStore, int count, int numQ, int *label, int *labelTable, long int *store, long int lookup, int type)
{
    int qa;
    int n = sqrt(numQ);
    double ran =((double)rand()/(double)RAND_MAX);
    int actual;
    int i;
    
    if (flag==0){
        qa = q-1;
        measurePerfectX(globalStore,q,numQ);
    
        if (type==1){
       
            int temp[10]={0};
            double ratio= (double)store[0]/(double)lookup;
            i=0;
            if (ratio<ran){
                
                while (ratio<ran){
                
                    i++;
                    ratio+=(double)store[i]/(double)lookup;
                }
                //
                int t=0;
                while(t<9) {
                    temp[t]=i%2;
                    i/= 2; t++;
                }
                //
                updateMeasure(temp,q,0,numQ,globalStore);
            }
        }
        
        actual = *(globalStore+2*qa+1);
    }
    
    if (flag==1){
        qa = q+1;
        measurePerfectZ(globalStore,q,numQ);
        
        if (type==1){
            
            int temp[10]={0};
            double ratio= (double)store[0]/(double)lookup;
            i=0;
            if (ratio<ran){
                
                while (ratio<ran){
                    
                    i++;
                    ratio+=(double)store[i]/(double)lookup;
                }
                //
                int t=0;
                while(t<9) {
                    temp[t]=i%2;
                    i/= 2; t++;
                }
                //
                updateMeasure(temp,q,1,numQ,globalStore);
            }
        }
        
        actual = *(globalStore+2*qa);
    }
    
    
    if (*(copy+2*qa)!=actual){
        
        *(labelTable+*label*3)=count;
        if (flag==0){
            *(labelTable+*label*3+1)=(qa-n)/(n*2);
            *(labelTable+*label*3+2)=qa%n/2;
        }
        if (flag==1){
            *(labelTable+*label*3+1)=(qa-1)/(n*2);
            *(labelTable+*label*3+2)=(qa-1)%n/2;
        }
        *label+=1;
        if (*label>1000) {printf("label beyond array size\n"); exit(1);}
    }
    *(copy+2*qa)=actual;
    //    sortError(globalStore,numQ);
}


void creatTable (int *labelTable, int (*table)[1000], int weightTime, int weightSpace, int numQ)
{
    int n=sqrt(numQ);
    int i=0; int j=0;
    
    int distanceTime; int distanceSpace; int distanceY; int distanceX;
    
    while (labelTable[i]!=0){
        
        j=3+i;
        
        while(labelTable[j]!=0){
            
            distanceTime=weightTime*(labelTable[j]-labelTable[i]);
            if ((labelTable[j+1]-labelTable[i+1])>=0){
                if ((labelTable[j+1]-labelTable[i+1])> (n/4))
                    distanceX = labelTable[i+1]+n/2-labelTable[j+1];
                else
                    distanceX = labelTable[j+1]-labelTable[i+1];
            }
            else{
                if ((labelTable[i+1]-labelTable[j+1])> (n/4))
                    distanceX = labelTable[j+1]+n/2-labelTable[i+1];
                else
                    distanceX = labelTable[i+1]-labelTable[j+1];
            }
            
            if ((labelTable[j+2]-labelTable[i+2])>=0){
                if ((labelTable[j+2]-labelTable[i+2])> (n/4))
                    distanceY = labelTable[i+2]+n/2-labelTable[j+2];
                else
                    distanceY = labelTable[j+2]-labelTable[i+2];
            }
            else{
                if ((labelTable[i+2]-labelTable[j+2])> (n/4))
                    distanceY = labelTable[j+2]+n/2-labelTable[i+2];
                else
                    distanceY= labelTable[i+2]-labelTable[j+2];
            }
            
            distanceSpace=weightSpace*(distanceX+distanceY);
            
            *(*(table+i/3)+j/3)=distanceTime+distanceSpace;
                       
            j+=3;
        }
        i+=3;
    }
}

int findLarger(int a, int b)
{
    if (a>b)
        return a;
    else
        return b;
}


int sortLogicalZ(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    int temp2=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[2*p+1];
    
    for (p=0;p<n;p+=2)
        temp2+=globalStore[2*(2*n-1+n*p)+1];
    
    flag=(temp1+temp2)%2;
    
    return flag;
    
}
void correctionX(int i, int j, int n, int *labelTableX, int *globalStore){

    {
        int p;
        int irow=labelTableX[3*i+1]; int icoll=labelTableX[3*i+2];
        int jrow=labelTableX[3*j+1]; int jcoll=labelTableX[3*j+2];
        
        if (irow!=jrow || icoll!=jcoll){
            
            
            int rowL = findLarger(irow,jrow);
            
            if (rowL!=jrow){
                int temprow=irow;
                int tempcoll=icoll;
                irow=jrow; icoll=jcoll;
                jrow=temprow; jcoll=tempcoll;
            }
            
            if ((jrow-irow)> n/4){
                for (p=0;p<(irow+n/2-jrow);p++){
                    int axis=irow*2*n+n+2*icoll-n*(2*p+1);
                    if (axis<0) axis+=numQ;
                    globalStore[2*axis+1]++;
                }
            }
            else{
                for (p=0;p<(jrow-irow);p++)
                    globalStore[2*(irow*2*n+n+2*icoll+n*(2*p+1))+1]++;
            }
            
            if (icoll>jcoll){
                if ((icoll-jcoll)> n/4){
                    for (p=0;p<(jcoll+n/2-icoll);p++){
                        int axis=jrow*2*n+n+2*jcoll-2*p-1;
                        if (axis<(jrow*2*n+n)) axis+=n;
                        globalStore[2*axis+1]++;
                    }
                }
                else{
                    for (p=0;p<(icoll-jcoll);p++)
                        globalStore[2*(jrow*2*n+n+2*jcoll+2*p+1)+1]++;
                }
            }
            if (icoll<jcoll){
                if ((jcoll-icoll)> n/4){
                    for (p=0;p<(icoll+n/2-jcoll);p++){
                        int axis=jrow*2*n+n+2*jcoll+2*p+1;
                        if (axis>(jrow*2*n+2*n-1)) axis-=n;
                        globalStore[2*axis+1]++;
                    }
                }
                else{
                    for (p=0;p<(jcoll-icoll);p++)
                        globalStore[2*(jrow*2*n+n+2*jcoll-2*p-1)+1]++;
                }
            }
            
        }
        
    }

}

void correctionZ(int i, int j, int n, int *labelTableZ, int *globalStore){
    
    int p;
    int irow=labelTableZ[3*i+1]; int icoll=labelTableZ[3*i+2];
    int jrow=labelTableZ[3*j+1]; int jcoll=labelTableZ[3*j+2];
    
    if (irow!=jrow || icoll!=jcoll){
        
        
        int rowL = findLarger(irow,jrow);
        
        if (rowL!=jrow){
            int temprow=irow;
            int tempcoll=icoll;
            irow=jrow; icoll=jcoll;
            jrow=temprow; jcoll=tempcoll;
        }
        //printf("irow=%d,icoll=%d,jrow=%d,jcoll=%d\n",irow,icoll,jrow,jcoll);
        if ((jrow-irow)> n/4){
            for (p=0;p<(irow+n/2-jrow);p++){
                int axis=irow*2*n+1+2*icoll-n*(2*p+1);
                if (axis<0) axis+=numQ;
                    globalStore[2*axis]++;
            }
        }
        else{
            for (p=0;p<(jrow-irow);p++)
                globalStore[2*(irow*2*n+1+2*icoll+n*(2*p+1))]++;
        }
        
        if (icoll>jcoll){
            if ((icoll-jcoll)> n/4){
                for (p=0;p<(jcoll+n/2-icoll);p++){
                    int axis=jrow*2*n+1+2*jcoll-2*p-1;
                    if (axis<(jrow*2*n)) axis+=n;
                        globalStore[2*axis]++;
                }
            }
            else{
                for (p=0;p<(icoll-jcoll);p++)
                    globalStore[2*(jrow*2*n+1+2*jcoll+2*p+1)]++;
            }
        }
        if (icoll<jcoll){
            if ((jcoll-icoll)> n/4){
                for (p=0;p<(icoll+n/2-jcoll);p++){
                    int axis=jrow*2*n+1+2*jcoll+2*p+1;
                    if (axis>(jrow*2*n+n-2)) axis-=n;
                        globalStore[2*axis]++;
                }
            }
            else{
                for (p=0;p<(jcoll-icoll);p++)
                    globalStore[2*(jrow*2*n+1+2*jcoll-2*p-1)]++;
            }
        }
        
    }
    
}

int sortLogicalX(int *globalStore, int numQ)
{
    
    int flag=0;
    int temp1=0;
    int temp2=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        temp1+=globalStore[2*(n*p)];
    
    for (p=0;p<n;p+=2)
        temp2+=globalStore[2*(n*(n-1)+p+1)];
    
    flag=(temp1+temp2)%2;
    
    return flag;
    
}
#endif
