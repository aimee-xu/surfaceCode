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

void cnot_v2 (int q1,int q2,int k,double *CnotPro,int *Q, int type){
    
    sortCnotError(q1,q2,Q);
    double ran;
    ran =((double)rand()/(double)RAND_MAX);
    
    if (type==1){
        
        double temp[7]={0};
        int i=0;
        temp[0]=CnotPro[8*k]+CnotPro[8*k+1];
        if (temp[0]<0.8)
            return;

        for (i=1;i<7;i++)
            temp[i]=temp[i-1]+CnotPro[8*k+i+1];

        if (ran>CnotPro[8*k] && ran<temp[0]){
            *(Q+2*q1+1)+=1;
        }
        else if (ran>temp[0] && ran<temp[1]){
            *(Q+2*q2)+=1;
        }
        else if (ran>temp[1] && ran<temp[2]){
            *(Q+2*q2)+=1;
            *(Q+2*q1+1)+=1;
        }
        else if (ran>temp[2] && ran<temp[3]){
            *(Q+2*q2+1)+=1;
            *(Q+2*q2)+=1;
        }
        else if (ran>temp[3] && ran<temp[4]){
            *(Q+2*q1+1)+=1;
            *(Q+2*q2+1)+=1;
            *(Q+2*q2)+=1;
        }
        else if (ran>temp[4] && ran<temp[5]){
            *(Q+2*q2+1)+=1;
        }
        else if (ran>temp[5] && ran<temp[6]){
            *(Q+2*q2+1)+=1;
        }
        else if (ran>temp[6] && ran<1){
            *(Q+2*q2+1)+=1; *(Q+2*q1+1)+=1;
        }
    
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

int findDeadQubit(int n, double *CnotPro, int *store, int *flagStore){

    int restart=0;
    int index; int p=0; int q=0;
    for (index=0;index<n/2;index++){
        int t1=2*n*index;
        int t2=index*2;
        int t3=2*n*(n/2-1)+2*index;
        int t4=(n/2-1)*2+2*n*index;
        if ((CnotPro[8*store[t1]]<0.8 || CnotPro[8*store[t2]]<0.8) || (CnotPro[8*store[t3]]<0.8 || CnotPro[8*store[t4]]<0.8)){
            restart=1;
            break;
        }
    }
    if (restart==1)
        return 1;
    
    for (p=0;p<(n/2)-1;p++){
        for (q=0;q<n/2-1;q++){
            int ts=n+1+2*n*p+2*q;
            if (CnotPro[8*store[ts]]<0.8){
                flagStore[4*(ts-n)+1]=ts;
                flagStore[4*(ts+n)]=ts;
                flagStore[4*(ts-1)+3]=ts;
                flagStore[4*(ts+1)+2]=ts;
            }
        }
    }
    
    for (p=1;p<(n/2)-1;p++){
        for (q=1;q<n/2-1;q++){
            int ts=2*n*p+2*q;
            if (CnotPro[8*store[ts]]<0.8){
                flagStore[4*(ts-n)+1]=ts;
                flagStore[4*(ts+n)]=ts;
                flagStore[4*(ts-1)+3]=ts;
                flagStore[4*(ts+1)+2]=ts;
            }
        }
    }

    return restart;
}


void parityCheckZ(double gateError, int *globalStore, int q3, int *store, double *CnotPro, int numQ, int direction, int type){
    
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
    
    measureNoise(qa,globalStore,gateError*type);
    cnot_v2 (q3,qa,store[q3],CnotPro,globalStore,type);
    cnot_v2 (q4,qa,store[q4],CnotPro,globalStore,type);
    if (direction!=0)
        cnot_v2 (q1,qa,store[q1],CnotPro,globalStore,type);
    if (direction!=1)
        cnot_v2 (q2,qa,store[q2],CnotPro,globalStore,type);
    measureNoise(qa,globalStore,gateError*type);
    
}

void parityCheckX(double gateError, int *globalStore, int q4, int *store, double *CnotPro, int numQ, int direction, int type){
    
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
    
    measureNoise(qa,globalStore,gateError*type);
    hadamard_r (qa,globalStore,0);
    cnot_v2 (qa,q1,store[q1],CnotPro,globalStore,type);
    cnot_v2 (qa,q2,store[q2],CnotPro,globalStore,type);
    if (direction!=0)
        cnot_v2 (qa,q3,store[q3],CnotPro,globalStore,type);
    if (direction!=1)
        cnot_v2 (qa,q4,store[q4],CnotPro,globalStore,type);
    
    hadamard_r (qa,globalStore,0);
    measureNoise(qa,globalStore,gateError*type);
    
}

void stabiliser(int q, int flag, int *copy, int *globalStore, int count, int numQ, int *label, int *labelTable, int *store, double *CnotPro,int type, int direction, double gateError)
{
    int qa;
    int n = sqrt(numQ);
    double ran =((double)rand()/(double)RAND_MAX);
    int actual;
    int i;
    
    if (flag==0){
        qa = q-1;
        parityCheckX(gateError,globalStore,q,store,CnotPro,numQ,direction,type);
    }
    
    if (flag==1){
        qa = q+1;
        parityCheckZ(gateError,globalStore,q,store,CnotPro,numQ,direction,type);
    }
    actual = *(globalStore+2*qa);
    if (*(copy+2*qa)!=actual){
        *(labelTable+*label*2)=count;
        *(labelTable+*label*2+1)=qa;
        *label+=1;
        if (*label>2000) {printf("label beyond array size\n"); exit(1);}
    }
    *(copy+2*qa)=actual;
}

int creatTable (int *labelTable, int label, int *table, int weightTime, int weightSpace, int numQ, int flag, int *flagStore)
{
    int n=sqrt(numQ);
    int i=0; int j=0;
    
    int distanceTime; int distanceSpace; int distanceY; int distanceX;
    int numError=0;
    while (i/2<label){
        
        int q1=labelTable[i+1];
        int q1x; int q1y; int q2x; int q2y;
        if (flag==0){
            q1x=(q1-n)/(n*2);
            q1y=q1%n/2;
        }
        if (flag==1){
            q1x=(q1-1)/(n*2);
            q1y=(q1-1)%n/2;
        }

        //to find out the mirror of the ith ancilla
        int distance;
        int temp;
        if (flag==0){
            distance=(n/2-1.5-q1x)*2;
            temp=(q1x+0.5)*2;
        }
        if (flag==1){
            distance=(n/2-1.5-q1y)*2;
            temp=(q1y+0.5)*2;
        }
        if (distance > temp)
            distance=temp;
        *(table+3*numError)=i/2; *(table+3*numError+1)=i/2+label;
        *(table+3*numError+2)=distance*weightSpace;
        numError++;
        
        j=2+i;
        while(j/2<label){//to find out every two pairs
            
            int q2=labelTable[j+1];
            
            if (flag==0){
                q2x=(q2-n)/(n*2);
                q2y=q2%n/2;
            }
            if (flag==1){
                q2x=(q2-1)/(n*2);
                q2y=(q2-1)%n/2;
            }

            int j1=((flagStore[4*q2]==flagStore[4*q1+1]) && (flagStore[4*q2]>0));
            int j2=((flagStore[4*q2+1]==flagStore[4*q1]) && (flagStore[4*q2+1]>0));
            int j3=((flagStore[4*q2+3]==flagStore[4*q1+2]) && (flagStore[4*q2+3]>0));
            int j4=((flagStore[4*q2+2]==flagStore[4*q1+3]) && (flagStore[4*q2+2]>0));
            
            if (j1 || j2 || j3 || j4){
                *(table+3*numError)=i/2; *(table+3*numError+1)=j/2;
                *(table+3*numError+2)=0;
            }
            else{
                
                distanceTime=weightTime*(labelTable[j]-labelTable[i]);
                
                if ((q2x-q1x)>=0)
                    distanceX = q2x-q1x;
                else
                    distanceX = q1x-q2x;
                
                if ((q2y-q1y)>=0)
                    distanceY = q2y-q1y;
                else
                    distanceY= q1y-q2y;
                
                distanceSpace=weightSpace*(distanceX+distanceY);
    
                *(table+3*numError)=i/2; *(table+3*numError+1)=j/2;
                *(table+3*numError+2)=distanceTime+distanceSpace;
            }
            numError++;
            j+=2;
        }
        i+=2;
    }
    for (i=label;i<2*label;i++){//to find out pairs between mirrors
        for (j=i+1;j<2*label;j++){
            *(table+3*numError)=i; *(table+3*numError+1)=j;
            *(table+3*numError+2)=0;
            numError++;
            
        }
    }
    return numError;
}

int findLarger(int a, int b)
{
    if (a>b)
        return a;
    else
        return b;
}

int findSmaller(int a, int b)
{
    if (a<b)
        return a;
    else
        return b;
}

void errorCorrect(int n,int jndex,int index,int label,int *labelTable, int type, int *globalStore){
    int p=0;
    
    if (type==0){
        
        if (jndex>=label && index>=label)
            ;
   
        else if (jndex>=label || index>=label){
            
            int smaller=findSmaller(jndex,index);
            
            int q1=labelTable[2*smaller+1];
            int row=(q1-n)/(n*2);
            int coll=q1%n/2;
            
            if (row<n/4){
                int axis=2*row*n+2*coll;
                while (axis>=0){
                    globalStore[2*axis+1]++;
                    axis=axis-2*n;
                }
            }
            else{
                int axis=(2*row+2)*n+2*coll;
                while (axis<=n*(n-2)+2*coll){
                    globalStore[2*axis+1]++;
                    axis=axis+2*n;
                }
            }
        }
        else{
            
            int q1=labelTable[2*index+1];
            int q2=labelTable[2*jndex+1];
            int irow=(q1-n)/(n*2); int icoll=q1%n/2;
            int jrow=(q2-n)/(n*2); int jcoll=q2%n/2;
            
            if (irow!=jrow || icoll!=jcoll){
                
                int rowL = findLarger(irow,jrow);
                
                if (rowL!=jrow){
                    int temprow=irow;
                    int tempcoll=icoll;
                    irow=jrow; icoll=jcoll;
                    jrow=temprow; jcoll=tempcoll;
                }
                
                for (p=0;p<(jrow-irow);p++)
                    globalStore[2*(irow*2*n+n+2*icoll+n*(2*p+1))+1]++;
                
                if (icoll>jcoll){
                    for (p=0;p<(icoll-jcoll);p++)
                        globalStore[2*(jrow*2*n+n+2*jcoll+2*p+1)+1]++;
                }
                
                if (icoll<jcoll){
                    for (p=0;p<(jcoll-icoll);p++)
                        globalStore[2*(jrow*2*n+n+2*jcoll-2*p-1)+1]++;
                }
            }
        }
    }
    if (type==1){
        
        if (jndex>=label && index>=label)
            ;
        else if (jndex>=label || index>=label){
            int smaller=findSmaller(jndex,index);
            int q1=labelTable[2*smaller+1];
            int row=(q1-1)/(n*2);
            int coll=(q1-1)%n/2;
            
            if (coll<n/4){
                int axis=2*n*row+coll*2;
                while (axis>=2*n*row){
                    globalStore[2*axis]++;
                    axis=axis-2;
                }
            }
            else{
                int axis=2*n*row+coll*2+2;
                while (axis<(2*row+1)*n){
                    globalStore[2*axis]++;
                    axis=axis+2;
                }
            }
        }
        else{
            int q1=labelTable[2*index+1];
            int q2=labelTable[2*jndex+1];
            int irow=(q1-1)/(n*2); int icoll=(q1-1)%n/2;
            int jrow=(q2-1)/(n*2); int jcoll=(q2-1)%n/2;
            if (irow!=jrow || icoll!=jcoll){
                
                int rowL = findLarger(irow,jrow);
                
                if (rowL!=jrow){
                    int temprow=irow;
                    int tempcoll=icoll;
                    irow=jrow; icoll=jcoll;
                    jrow=temprow; jcoll=tempcoll;
                }
                
                for (p=0;p<(jrow-irow);p++)
                    globalStore[2*(irow*2*n+1+2*icoll+n*(2*p+1))]++;
                
                if (icoll>jcoll){
                    for (p=0;p<(icoll-jcoll);p++)
                        globalStore[2*(jrow*2*n+1+2*jcoll+2*p+1)]++;
                }
                if (icoll<jcoll){
                    for (p=0;p<(jcoll-icoll);p++)
                        globalStore[2*(jrow*2*n+1+2*jcoll-2*p-1)]++;
                }
            }
            
        }
        
    }
}

int sortLogicalZ(int *globalStore, int numQ)
{
    
    int flag=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        flag+=globalStore[2*p+1];

    flag=flag%2;
    
    return flag;
    
}

int sortLogicalX(int *globalStore, int numQ)
{
    
    int flag=0;
    
    int p; int q; int n =sqrt(numQ);
    
    for (p=0;p<n;p+=2)
        flag+=globalStore[2*(n*p)];

    flag=flag%2;
    
    return flag;
    
}
