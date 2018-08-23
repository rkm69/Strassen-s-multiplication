/*####################  MATRIX MULTIPLICATION USING STRASSEN's METHOD   #########################*/
        //Author:RAHUL KUMAR
        //CSE student
        //TECHNO INDIA,Saltlake


#include<stdio.h>
#include<math.h>
typedef struct {
    int mat[16][16];
    int size;
}matrix;

typedef struct {
    int row,col;
}row_col;

row_col read_rowcol(){
    row_col r1;
    int r,c;
    printf("Enter the row and column::");
    scanf("%d %d",&r,&c);
    r1.row=r;
    r1.col=c;
    return r1;
}


matrix read(row_col p){
    int i,j;
    matrix mat1;
    printf("Enter the elements of matrix::");
    for(i=0;i<p.row;i++){
        for(j=0;j<p.col;j++){
            scanf("%d",&(mat1.mat[i][j]));
        }
    }
    mat1.size=2;
    if(p.col>2||p.row>2)
        mat1.size=4;
    else if(p.col>4||p.row>4)
        mat1.size=8;
    else if(p.col>8||p.row>8)
        mat1.size=16;
    for(i=0;i<p.row;i++){
        for(j=p.col;j<mat1.size;j++){
            mat1.mat[i][j]=0;
        }
    }
    for(i=p.row;i<mat1.size;i++){
        for(j=0;j<mat1.size;j++){
            mat1.mat[i][j]=0;
        }
    }
    
    return mat1;
    
}

void print(matrix mat1){
    int i,j;
    for(i=0;i<mat1.size;i++){
        for(j=0;j<mat1.size;j++){
            printf("%5d ",mat1.mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

matrix add(matrix m1,matrix m2){
    matrix m3;
    int i,j;
    m3.size=m1.size;
    for(i=0;i<m1.size;i++){
        for(j=0;j<m1.size;j++){
            m3.mat[i][j]=m1.mat[i][j]+m2.mat[i][j];
        }
    }
    return m3;
}

matrix sub(matrix m1,matrix m2){
    matrix m3;
    int i,j;
    m3.size=m1.size;
    for(i=0;i<m1.size;i++){
        for(j=0;j<m1.size;j++){
            m3.mat[i][j]=m1.mat[i][j]-m2.mat[i][j];
        }
    }
    return m3;
}


matrix strassen(matrix a,matrix b,int k){
    
    if(k==2){
        int P,Q,R,S,T,U,V;
        P=(a.mat[0][0]+a.mat[1][1])*(b.mat[0][0]+b.mat[1][1]);
        Q=(a.mat[1][0]+a.mat[1][1])*b.mat[0][0];
        R=a.mat[0][0]*(b.mat[0][1]-b.mat[1][1]);
        S=a.mat[1][1]*(b.mat[1][0]-b.mat[0][0]);
        T=(a.mat[0][0]+a.mat[0][1])*b.mat[1][1];
        U=(a.mat[1][0]-a.mat[0][0])*(b.mat[0][0]+b.mat[0][1]);
        V=(a.mat[0][1]-a.mat[1][1])*(b.mat[1][0]+b.mat[1][1]);
        matrix r;
        r.size=2;
        r.mat[0][0]=P+S-T+V;
        r.mat[0][1]=R+T;
        r.mat[1][0]=Q+S;
        r.mat[1][1]=P-Q+R+U;
        return r;

    }
    else{
    k=k/2;
    int i,j;
    matrix A11,A12,A21,A22,  B11,B12,B21,B22,  p,q,r,s,t,u,v,   C11,C12,C21,C22;
    A11.size=A12.size=A21.size=A22.size=B11.size=B12.size=B21.size=B22.size=k;
    for(i=0;i<k;i++){
        for(j=0;j<k;j++){
            A11.mat[i][j]=a.mat[i][j];
            B11.mat[i][j]=b.mat[i][j];
            A12.mat[i][j]=a.mat[i][k+j];
            B12.mat[i][j]=b.mat[i][k+j];
            A21.mat[i][j]=a.mat[k+i][j];
            B21.mat[i][j]=b.mat[k+i][j];
            A22.mat[i][j]=a.mat[k+i][k+j];
            B22.mat[i][j]=b.mat[k+i][k+j];
        }
    }
    printf("The partion of the Initial Matrix1, whose size is:: %d\n",k); 
    printf("A11::\n");
    print(A11);
    printf("A12::\n");    
    print(A12);
    printf("A21::\n");
    print(A21);
    printf("A22::\n");
    print(A22);
    printf("The partion of the Initial Matrix1, whose size is:: %d\n",k); 
    printf("B11::\n");
    print(B11);
    printf("B12::\n");
    print(B12);
    printf("B21::\n");
    print(B21);
    printf("B22::\n");
    print(B22);
    
    printf("The values of the 7 constant matrices p,q,r,s,t,u,v are as follows each of size %d\n",k);
    
    p=strassen(add(A11,A22),add(B11,B22),k);  //p
    printf("p::\n");
    print(p); 
    q=strassen(add(A21,A22),B11,k);  //q
    printf("q::\n");
    print(q);
    r=strassen(A11,sub(B12,B22),k); //r
    printf("r::\n");
    print(r);
    s=strassen(A22,sub(B21,B11),k); //s
    printf("s::\n");
    print(s);
    t=strassen(add(A11,A12),B22,k); //t
    printf("t::\n");
    print(t);
    u=strassen(sub(A21,A11),add(B11,B12),k);    //u
    printf("u::\n");
    print(u);
    v=strassen(sub(A12,A22),add(B21,B22),k);  //v      
    printf("v::\n");
    print(v);
    printf("The values of the partioned value of resultant matrix each of size %d\n",k);
    C11=add(sub(p,t),add(s,v)); 
    printf("C11::\n");
    print(C11);
    C12=add(r,t);
    printf("C12::\n");
    print(C12);
    C21=add(q,s);
    printf("C21::\n");
    print(C21);
    C22=add(sub(p,q),add(r,u));
    printf("C22::\n");
    print(C22);
    matrix res;
    res.size=k*2;
    for(i=0;i<k;i++){
        for(j=0;j<k;j++){
            res.mat[i][j]=C11.mat[i][j];
            res.mat[i][k+j]=C12.mat[i][j];
            res.mat[k+i][j]=C21.mat[i][j];
            res.mat[k+i][k+j]=C22.mat[i][j];
        }
    }
    return res;
        
    }
   
    
}


int main(){
    matrix m1,m2,m3;
    row_col r1,r2;
    printf("MAtrix 1::");
    r1=read_rowcol();
    
    m1=read(r1);
    printf("matrix1::\n");
    print(m1);
    printf("MAtrix 2::");
    r2=read_rowcol();
    
    m2=read(r2);
    printf("Matrix2::\n");
    print(m2);
    m3=strassen(m1,m2,m1.size);
    printf("The multiplication using strassen's method\n");
    printf("The resultant MATRIX C[%d][%d] is ::\n",m1.size,m1.size);
    print(m3);
    
}
    


    

