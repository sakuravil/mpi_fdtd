#define X_m_max 3200
#define Y_m_max 3200
#define T_max_m 5000
#define MAX_STRING 128
#define PI 3.14159265358979

static  double  P1[X_m_max][Y_m_max];
static  double  Vx[X_m_max][Y_m_max];
static  double  Vy[X_m_max][Y_m_max];
static  double  send[3200];
static  double  recv[3200]={0},recv_P1[X_m_max][Y_m_max];
static  double  data_s[3200]={0};
static  double  data_r[3200]={0};
static  double  PX1[6][Y_m_max];
static  double  PY1[6][X_m_max];
static  double  PX2[6][Y_m_max];
static  double  PY2[6][X_m_max];
static  double  amp1[T_max_m],amp2[T_max_m];
static  double  time_e[10000]={0},time_s[10000]={0};

static  double  f, dx, dt, range, depth, cal_time;
static  double  rou0, c0, gensui0, absp0, alpha0, kap0, hasu0, c_m0, hasu_o0;
static  double  Cp1, Cp2, Cv1, Cv2;
static  double  Ca0, Ca1, Ca2, Cah1, Cah2, a1, a2, d1, d2;

static  int     tag=0,id=0,start_0,start,fin,idd,pro,rank,name_len,i, j, T, del_T, l, m, n, x_max, y_max, T_max ,k;
static  int     i_min, i_max, j_min, j_max, W_end, WN , end;
static  int     ng,mg;
static  int     SX, SY, RX, RY, SC, BC;
static  double  sd, rd, bd, rds;
static  int     wide_S,high_S,wide_E,high_E,wd,hg;
static  double  stopwatch[3][10]={};
static  int     Process_Point_X,Process_Point_Y; 
static  int     P;
static  double  set[10000]={};

void gatherdata(int);

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>
#include        <omp.h>
#include        "fdtd_func.h"
#include        "boundary.h"
#include        "mpi.h"

void Split(int x,int y,int id){
  int num[10000]={};
  int count=0;

  // prime factorization
  for(i=1;i<=P;i++)
    if(P%i==0)num[count++]=i;

  if(count%2==0){
    wd=num[count/2-1];
    hg=num[count/2];
  } else{
    wd=num[count/2];
    hg=num[count/2];
  }
  // Calculate calculation range
  high_S=y/hg*(id%hg);
  high_E=y/hg*(id%hg+1);
  wide_S=x/wd*(id/hg);
  wide_E=x/wd*(id/hg+1);

  /* Display placement coordinates */
  Process_Point_X=wd*wide_S/x_max;
  Process_Point_Y=hg*high_S/y_max;
  printf("[P %d]PROCESS(x,y)=(%d,%d)\n",id,Process_Point_X,Process_Point_Y);

  if(Process_Point_X!=0)wide_S--;
  if(Process_Point_Y!=0)high_S--;
  wide_E--;
  high_E--;

  printf("high_S:%d high_E:%d\n",high_S,high_E);
  printf("wide_S:%d wide_E:%d\n",wide_S,wide_E);
}

void Pressure(int id){
  if( wide_S<=SX && wide_E>=SX && high_S<=SY && high_E>=SY && T <= W_end)P1[SX][SY]=Driver(T);
  #pragma omp parallel for private(j)
  for(i = 1 + wide_S; i <= wide_E; i+=2){
    for(j = 1 + high_S; j <= high_E; j+=2){
      P1[i][j] = Cp1*P1[i][j] - Cp2*(Vx[i][j]-Vx[i-1][j]+Vy[i][j]-Vy[i][j-1]);
      P1[i][j+1] = Cp1*P1[i][j+1] - Cp2*(Vx[i][j+1]-Vx[i-1][j+1]+Vy[i][j+1]-Vy[i][j]);

      P1[i+1][j] = Cp1*P1[i+1][j] - Cp2*(Vx[i+1][j]-Vx[i][j]+Vy[i+1][j]-Vy[i+1][j-1]);
      P1[i+1][j+1] = Cp1*P1[i+1][j+1] - Cp2*(Vx[i+1][j+1]-Vx[i][j+1]+Vy[i+1][j+1]-Vy[i+1][j]);
    }
  }
}

void Velocity(int id){
  #pragma omp parallel
  {
    #pragma omp for private(j)
    for(i = wide_S; i <= wide_E; i+=2){
      for(j = 1 + high_S; j <= high_E; j++){
        Vx[i][j] = Cv1*Vx[i][j] - Cv2*(P1[i+1][j]-P1[i][j]);
        Vx[i+1][j] = Cv1*Vx[i+1][j] - Cv2*(P1[i+2][j]-P1[i+1][j]);
      }
    }

    #pragma omp for private(j)
    for(i = 1 + wide_S; i <= wide_E; i++){
      for(j = high_S ; j <= high_E ; j+=2){
        Vy[i][j] = Cv1*Vy[i][j] - Cv2*(P1[i][j+1]-P1[i][j]);
        Vy[i][j+1] = Cv1*Vy[i][j+1] - Cv2*(P1[i][j+2]-P1[i][j+1]);
      }
    }
  }
}

void send_setting(MPI_Request send_req[10000],MPI_Request recv_req[10000]){
  // Not located on the left end
  if(Process_Point_X != 0){
    MPI_Recv_init(data_r,y_max/hg,MPI_DOUBLE,id-hg,0,MPI_COMM_WORLD,&recv_req[id-hg]);
    MPI_Send_init(data_s,y_max/hg,MPI_DOUBLE,id-hg,1,MPI_COMM_WORLD,&send_req[id-hg]);
  }
  // Not located on the right end
  if(Process_Point_X != wd-1){
    MPI_Send_init(data_s,y_max/hg,MPI_DOUBLE,id+hg,0,MPI_COMM_WORLD,&send_req[id+hg]);
    MPI_Recv_init(data_r,y_max/hg,MPI_DOUBLE,id+hg,1,MPI_COMM_WORLD,&recv_req[id+hg]);
  }

  // Not located on the upper end
  if(Process_Point_Y != 0){
    MPI_Recv_init(data_r,x_max/wd,MPI_DOUBLE,id-1,2,MPI_COMM_WORLD,&recv_req[id-1]);
    MPI_Send_init(data_s,x_max/wd,MPI_DOUBLE,id-1,3,MPI_COMM_WORLD,&send_req[id-1]);
  }

  // Not located on the lower end
  if(Process_Point_Y != hg-1){
    MPI_Send_init(data_s,x_max/wd,MPI_DOUBLE,id+1,2,MPI_COMM_WORLD,&send_req[id+1]);
    MPI_Recv_init(data_r,x_max/wd,MPI_DOUBLE,id+1,3,MPI_COMM_WORLD,&recv_req[id+1]);
  }
}

void share(int id,MPI_Status status,MPI_Request send_req[10000],MPI_Request recv_req[10000]){
  int count=0;
  double ss,ee;

  // Not located on the left end
  if(Process_Point_X != 0){
    ss=MPI_Wtime();
    MPI_Start(&recv_req[id-hg]);
    MPI_Wait(&recv_req[id-hg],&status); 
    ee=MPI_Wtime();
    set[1]+=ee-ss;
    for(i = 1 + high_S; i <= high_E; ++i,++count) P1[wide_S][i]=data_r[count];
    count=0;

    for(i = 1 + high_S; i <= high_E; ++i,++count) data_s[count]=P1[wide_S+1][i]; 
    ss=MPI_Wtime();
    MPI_Start(&send_req[id-hg]);
    MPI_Wait(&send_req[id-hg],&status);
    ee=MPI_Wtime();
    set[0]+=ee-ss;

    count=0;
  }

  // Not located on the right end
  if( Process_Point_X != wd-1){
    for(i = 1 + high_S; i <= high_E; ++i,++count) data_s[count]=P1[wide_E][i];
    ss=MPI_Wtime(); 
    MPI_Start(&send_req[id+hg]);
    MPI_Wait(&send_req[id+hg],&status);
    ee=MPI_Wtime();
    set[0]+=ee-ss;
    count=0;

    ss=MPI_Wtime();
    MPI_Start(&recv_req[id+hg]);
    MPI_Wait(&recv_req[id+hg],&status);
    ee=MPI_Wtime();
    set[1]+=ee-ss;
    for(i = 1 + high_S; i <= high_E; ++i,++count) P1[wide_E+1][i]=data_r[count];
    count=0;
  }

  // Not located on the upper end
  if(Process_Point_Y != 0){
    ss=MPI_Wtime();
    MPI_Start(&recv_req[id-1]);
    MPI_Wait(&recv_req[id-1],&status);
    ee=MPI_Wtime();
    set[1]+=ee-ss;
    for(i = 1 + wide_S; i <= wide_E; ++i,++count) P1[i][high_S]=data_r[count];
    count=0;

    for(i = 1 + wide_S; i <= wide_E; ++i,++count) data_s[count]=P1[i][high_S+1]; 
    ss=MPI_Wtime();
    MPI_Start(&send_req[id-1]);
    MPI_Wait(&send_req[id-1],&status);
    ee=MPI_Wtime();
    set[0]+=ee-ss;
    count=0;
  }

  // Not located on the lower end
  if( Process_Point_Y != hg-1){
    for(i = 1 + wide_S; i <= wide_E; ++i,++count) data_s[count]=P1[i][high_E];
    ss=MPI_Wtime();
    MPI_Start(&send_req[id+1]);
    MPI_Wait(&send_req[id+1],&status);
    ee=MPI_Wtime();
    set[0]+=ee-ss;
    count=0;

    ss=MPI_Wtime();
    MPI_Start(&recv_req[id+1]);
    MPI_Wait(&recv_req[id+1],&status);
    ee=MPI_Wtime();
    set[1]+=ee-ss;
    for(i = 1 + wide_S; i <= wide_E; ++i,++count)P1[i][high_E+1]=data_r[count];
    count=0;
  }
}


  MPI_Status status;
  int count=0;

  for(i = 1 ; i <= x_max/wd ; i++){
    if(id!=P-1){
      for(j = 1 + high_S; j <= high_E; ++j,++count) data_s[count]=P1[i + wide_S][j];
      MPI_Send(data_s,y_max/hg,MPI_DOUBLE,P-1,0,MPI_COMM_WORLD);
      count=0;
    } else{
      for(j=0;j<P-1;j++){
        MPI_Recv(data_r,y_max/hg,MPI_DOUBLE,j,0,MPI_COMM_WORLD,&status);
        for(k = 1 + y_max/hg*(j%hg); k <= y_max/hg*(j%hg+1) ; ++k,++count)
          P1[i + x_max/wd*(j/hg)][k]=data_r[count];
        count=0;
      }
    }
  }
}

void snap(int id){
  gatherdata(id);
  if(rank==P-1) FieldSnap(P1);
}

void logdata(int id,int T){
  if(T%(T_max/10)==0 && P1[x_max/4][y_max/2]!=0)
    printf("1   %d:P1=%f  rank:%d\n",T/10,10.0*log10(P1[x_max/4][y_max/2]*P1[x_max/4][y_max/2]),id);
  if(T%(T_max/10)==0 && P1[x_max/2][y_max/2]!=0)    
    printf("2   %d:P1=%f  rank:%d\n",T/10,10.0*log10(P1[x_max/2][y_max/2]*P1[x_max/2][y_max/2]),id);
  if(T%(T_max/10)==0 && P1[x_max/4*3][y_max/2]!=0)
    printf("3   %d:P1=%f  rank:%d\n",T/10,10.0*log10(P1[x_max/4*3][y_max/2]*P1[x_max/4*3][y_max/2]),id);
}

int main(int argc,char** argv){   
  static char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;
  MPI_Request send_req[10000],recv_req[10000];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name,&name_len); 

  time_s[rank]=MPI_Wtime();

  for(i=0;i<P;i++)if(rank==i)id = i;

  Make_param();
  Initialize(id);
  Split(x_max,y_max,id);
  send_setting(send_req,recv_req);

  if(rank==0){
    printf("\n**************** 2D FD-TD *******************\n");
    printf("Number of Process:%d\n",P);
    printf("   Starting simulation of sound propagation\n\n");
    printf(" \n Model parameter\n");
    printf("   Frequency = %f[Hz]\n",(float)f);
    printf("   Range = %f[m]\n",(float)range);
    printf("   Depth = %f[m]\n",(float)depth);
    printf(" \n Grid parameter\n");
    printf("   dx = %f[mm]\n",(float)dx*1000);
    printf("   dt = %f[ms]\n\n",(float)dt*1000);
  }

  for(T = 0; T <= T_max; ++T){
    stopwatch[0][0]=MPI_Wtime();

    Pressure(id);
    stopwatch[0][1]=MPI_Wtime();
    boundary(id);
    stopwatch[0][2]=MPI_Wtime();
    share(id,status,send_req,recv_req);
    stopwatch[0][3]=MPI_Wtime();
    Velocity(id);
    stopwatch[0][4]=MPI_Wtime();
    WaveT(id);
    //if( T%del_T==0 && T >= del_T)snap(id);
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
    OutputF();///計算結果出力
    printf("\n  calculation complete !\n\n");
  }
  OutputW(id);///送波・受波出力

  time_e[rank]=MPI_Wtime();
  printf("rank:%d time:%f[s]\n",rank,time_e[rank]-time_s[rank]);///各プロセスのかかった時間

  if(id==P-1){
    printf("Pressure(): %f[s]\n",stopwatch[1][1]);
    printf("boundary(): %f[s]\n",stopwatch[1][2]);
    printf("share(): %f[s]\n",stopwatch[1][3]);
    printf("...[1]: %f[s]\n",set[0]);
    printf("...[2]: %f[s]\n",set[1]);
    printf("Velocity(): %f[s]\n",stopwatch[1][4]);
    printf("gatherdata(): %f[s]\n",stopwatch[1][5]);
  }
  MPI_Finalize(); 
}
