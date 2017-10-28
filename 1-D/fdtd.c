#define X_m_max 3200
#define Y_m_max 3200
#define T_max_m 5000
#define MAX_STRING 128
#define PI 3.14159265358979

static  float  P1[X_m_max][Y_m_max];
static  float  Vx[X_m_max][Y_m_max];
static  float  Vy[X_m_max][Y_m_max];
static  float  send[3200];
static  float  recv[3200]={0},recv_P1[X_m_max][Y_m_max];
static  float  data_s[3200]={0};
static  float  data_r[3200]={0};
static  float  PX1[6][Y_m_max];
static  float  PY1[6][X_m_max];
static  float  PX2[6][Y_m_max];
static  float  PY2[6][X_m_max];
static  float  amp1[T_max_m],amp2[T_max_m];
static  float  time_e[10000]={0},time_s[10000]={0};
static  float  stopwatch[3][10]={};

static  float  f, dx, dt, range, depth, cal_time;
static  float  rou0, c0, gensui0, absp0, alpha0, kap0, hasu0, c_m0, hasu_o0;
static  float  Cp1, Cp2, Cv1, Cv2;
static  float  Ca0, Ca1, Ca2, Cah1, Cah2, a1, a2, d1, d2;

static  int     tag=0,id=0,start_0,start,fin,idd,pro,rank,name_len,i, j, T, del_T, l, m, n, x_max, y_max, T_max ,k;
static  int     i_min, i_max, j_min, j_max, W_end, WN , end;
static  int     ng,mg;
static  int     SX, SY, RX, RY, SC, BC;
static  float  sd, rd, bd, rds;
static  int     P;
static  float  set[10000]={};
static  float  ss,ee;

void gatherdata(int);

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>
#include        <omp.h>
#include        "fdtd_func.h"
#include        "boundary.h"
#include        "mpi.h"

void Pressure(int id){
  start=((x_max)/P)*id;

  if( start<=SX && x_max/P + start>=SX && T <= W_end)P1[SX][SY]=Driver(T);
  #pragma omp parallel for private(j)
  for(i = 1 + start; i <= x_max/P + start; i+=2){
    for(j = 1; j <= y_max-2; j+=2){
      P1[i][j] = Cp1*P1[i][j] - Cp2*(Vx[i][j]-Vx[i-1][j]+Vy[i][j]-Vy[i][j-1]);
      P1[i][j+1] = Cp1*P1[i][j+1] - Cp2*(Vx[i][j+1]-Vx[i-1][j+1]+Vy[i][j+1]-Vy[i][j]);

      P1[i+1][j] = Cp1*P1[i+1][j] - Cp2*(Vx[i+1][j]-Vx[i][j]+Vy[i+1][j]-Vy[i+1][j-1]);
      P1[i+1][j+1] = Cp1*P1[i+1][j+1] - Cp2*(Vx[i+1][j+1]-Vx[i][j+1]+Vy[i+1][j+1]-Vy[i+1][j]);
    }
  }
}

void Velocity(int id){

  start=((x_max)/P)*id;
  #pragma omp parallel
  {
  #pragma omp for private(j)
  for(i = start; i <= x_max/P + start; i+=2){
    for(j = 1; j <= y_max - 2; j++){
      Vx[i][j] = Cv1*Vx[i][j] - Cv2*(P1[i+1][j]-P1[i][j]);
      Vx[i+1][j] = Cv1*Vx[i+1][j] - Cv2*(P1[i+2][j]-P1[i+1][j]);
    }
  }

  #pragma omp for private(j)
  for(i = 1 + start; i <= x_max/P + start + 1; i++){
    for(j = 0; j <= y_max - 2; j+=2){
      Vy[i][j] = Cv1*Vy[i][j] - Cv2*(P1[i][j+1]-P1[i][j]);
      Vy[i][j+1] = Cv1*Vy[i][j+1] - Cv2*(P1[i][j+2]-P1[i][j+1]);
    }
  }
  }
}

void send_setting(MPI_Request req[1000],MPI_Request rreq[1000]){
  if(id!=0 && id<P){
    MPI_Send_init(data_s,y_max-1,MPI_float,id-1,0,MPI_COMM_WORLD,&req[id-1]);
    MPI_Recv_init(data_r,y_max-1,MPI_float,id-1,0,MPI_COMM_WORLD,&rreq[id-1]);
  }

  if(id!=P-1 && id<P){
    MPI_Send_init(data_s,y_max-1,MPI_float,id+1,0,MPI_COMM_WORLD,&req[id+1]);
    MPI_Recv_init(data_r,y_max-1,MPI_float,id+1,0,MPI_COMM_WORLD,&rreq[id+1]);
  }
}

void share(int id,MPI_Status status,MPI_Request req[1000],MPI_Request rreq[1000]){
  if(id!=0 && id<P){
    ss=MPI_Wtime();
    MPI_Start(&rreq[id-1]);
    MPI_Wait(&rreq[id-1],&status);
    ee=MPI_Wtime();
    set[0]+=ee-ss;

    ss=MPI_Wtime();
    for(i = 1; i <= y_max-1; ++i)P1[x_max/P*id][i]=data_r[i];
    ee=MPI_Wtime();
    set[1]+=ee-ss;

    ss=MPI_Wtime();
    for(i = 1; i <= y_max-1; ++i)data_s[i]=P1[x_max/P*id+1][i]; 
    ee=MPI_Wtime();
    set[2]+=ee-ss;

    ss=MPI_Wtime();
    MPI_Start(&req[id-1]);
    MPI_Wait(&req[id-1],&status);
    ee=MPI_Wtime();
    set[3]+=ee-ss;
  }

  if(id!=P-1 && id<P){
    for(i = 1; i <= y_max-1; ++i)data_s[i]=P1[x_max/P*(id+1)][i]; 
    MPI_Start(&req[id+1]);
    MPI_Wait(&req[id+1],&status);
    MPI_Start(&rreq[id+1]);
    MPI_Wait(&rreq[id+1],&status); 
    for( i = 1; i <= y_max-1 ; ++i )P1[x_max/P*(id+1)+1][i]=data_r[i];
  }
}

void gatherdata(int id){
  MPI_Status status;

  for(i = 1 ;i <= (x_max)/P ;i++){
    if(id!=P-1){
      #pragma omp parallel for 
      for( j = 1; j <= y_max-1 ; ++j )data_s[j]=P1[i+x_max/P*id][j]; 
        MPI_Send(data_s,y_max-1,MPI_float,P-1,0,MPI_COMM_WORLD);      
    } else{
        for(j=0;j<P-1;j++){
          MPI_Recv(data_r,y_max-1,MPI_float,j,0,MPI_COMM_WORLD,&status);
          for( k = 1; k <= y_max-1 ; ++k )P1[i+x_max/P*j][k]=data_r[k];
        }
    }
  }
}

void snap(int id){
  gatherdata(id);
  if(rank==P-1)FieldSnap(P1);
}

void logdata(int id,int T){
  if(T%(T_max/10)==0 && P1[x_max/4][y_max/2]!=0)
    printf("1   %d:P1=%f\n",T/10,10.0*log10(P1[x_max/4][y_max/2]*P1[x_max/4][y_max/2]));
  if(T%(T_max/10)==0 && P1[x_max/2][y_max/2]!=0)    
    printf("2   %d:P1=%f\n",T/10,10.0*log10(P1[x_max/2][y_max/2]*P1[x_max/2][y_max/2]));
  if(T%(T_max/10)==0 && P1[x_max/4*3][y_max/2]!=0)
    printf("3   %d:P1=%f\n",T/10,10.0*log10(P1[x_max/4*3][y_max/2]*P1[x_max/4*3][y_max/2]));
}

int main(int argc,char** argv){   
  static char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;
  MPI_Request req[1000],rreq[1000];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name,&name_len); 
  time_s[rank]=MPI_Wtime();

  for(i=0;i<P;i++)if(rank==i)id = i;

  Make_param();
  Initialize(id);
  send_setting(req,rreq);

  if(rank==0){
    printf("\n**************** 2D FD-TD *******************\n");
    printf("   Starting simulation of sound propagation\n\n");
    printf(" \n Model parameter\n");
    printf("   Frequency = %f[Hz]\n",(float)f);
    printf("   Range = %f[m]\n",(float)range);
    printf("   Depth = %f[m]\n",(float)depth);
    printf(" \n Grid parameter\n");
    printf("   dx = %f[mm]\n",(float)dx*1000);
    printf("   dt = %f[ms]\n\n",(float)dt*1000);
  }
  omp_set_num_threads(2);
  
  // main calculation 
  for( T = 0 ; T <= T_max ; ++T ){
    stopwatch[0][0]=MPI_Wtime();
    Pressure(id);
    stopwatch[0][1]=MPI_Wtime();
    boundary(id);
    stopwatch[0][2]=MPI_Wtime();
    share(id,status,req,rreq);
    stopwatch[0][3]=MPI_Wtime();
    Velocity(id);
    stopwatch[0][4]=MPI_Wtime();
    WaveT(id);
    logdata(id,T);
    //if(T%(T_max/50) == 0)
    //printf("rank=%d  %4d%% calculation complete\n",rank,T*100/T_max);
    for(i=1;i<5;i++)stopwatch[1][i]+=(stopwatch[0][i]-stopwatch[0][i-1]);
  }

  stopwatch[0][0]=MPI_Wtime();
  gatherdata(id);
  stopwatch[0][5]=MPI_Wtime();
  stopwatch[1][5]=stopwatch[0][5]-stopwatch[0][0];

  if(rank==P-1){
    OutputF();
    printf("\n  calculation complete !\n\n");
  }
  OutputW(id);
  time_e[rank]=MPI_Wtime();
  printf("rank:%d time:%f[s]\n",rank,time_e[rank]-time_s[rank]);

  if(id==P-1){
    printf("Pressure(): %f[s]\n",stopwatch[1][1]);
    printf("boundary(): %f[s]\n",stopwatch[1][2]);
    printf("share(): %f[s]\n",stopwatch[1][3]);
    printf("...[1]: %f[s]\n",set[0]);
    printf("...[2]: %f[s]\n",set[1]);
    printf("...[3]: %f[s]\n",set[2]);
    printf("...[4]: %f[s]\n",set[3]);
    printf("Velocity(): %f[s]\n",stopwatch[1][4]);
    printf("gatherdata(): %f[s]\n",stopwatch[1][5]);
  }
  MPI_Finalize(); 
}
