#ifndef FDTD_FUNC_H
#define FDTD_FUNC_H

void Initialize(int id){
  for(i = wide_S; i <= wide_E; i++){
    for(j = high_S; j <= high_E; ++j){
      P1[i][j] = 0;
      Vx[i][j] = 0;
      Vy[i][j] = 0;
    }
  }
}

void Make_param(void){
  int snap;

  // parameter input 
  range = 1200;   // propagation distance [m] 
  depth = 1200;   // depth [m]
  cal_time = 1.2; // propagation time [sec]

  dx = 3.0;       // space step width [m]
  dt = 0.0006;    // time step width [sec]

  f = 25;         // excitation frequency
  WN = 4;         // excitation pulse length 

  // SC = 0: continuous wave,
  //      1: gauusian pulse wave,
  //      2: hanning pulse wave
  SC = 2;         // excitation function

  rou0 = 1000;    // density
  c0 = 1500;      // sound velocity
  absp0 = 0.0;    // absorption coefficient [dB/λ] 
  gensui0 = (f*absp0)/(8.686*c0); // damping coefficient [neper/m]

  // BC = 0: Mur_1st, <<<<<<<------only use
  //      1: Mur_2nd,
  //      2: Higdon_1st, 
  //      3: Higdon_2nd
  BC = 0;         // absorption boundary condition

  x_max = (int)(range/dx);
  y_max = (int)(depth/dx);
  T_max = (int)(cal_time/dt);

  SX = (int)(x_max/2);    // transmission point (X)
  SY = (int)(y_max/2);    // transmission point (Y)

  RX = (int)(4*x_max/5);  // receiving point (X)
  RY = (int)(y_max/4);    // receiving point (Y)

  snap = 10;              // number of snapshot  
  del_T = (int)(T_max/snap);


  if(SC== 0)W_end = T_max;  // float the excitation finish step count
  else if (SC == 1)W_end = (int)(2*WN/(f * dt));
  else if (SC == 2)W_end = (int)(WN/(f * dt));

  // file thinning count
  mg = 0;
  ng = 0;

  // wave number including loss,
  // sound velocity, 
  // cofficient of loss term      caluculation
  hasu_o0 = 2*PI*f/c0;                                  // actual wave number 
  hasu0 = sqrt(hasu_o0*hasu_o0 - gensui0*gensui0);      // wave number in case of loss
  c_m0 = 2*PI*f/hasu0;          
  alpha0 = 2*hasu_o0*gensui0*rou0*c_m0/hasu0;           // absorption term η

  // caluculation of difference expression coefficient
  kap0 = c_m0 * c_m0 * rou0;
  Cp1 = 1.0;
  Cp2 = kap0 * dt / dx;
  Cv1 = (2*rou0 - alpha0*dt)/(2*rou0 + alpha0*dt);
  Cv2 = dt*2/((2*rou0 + alpha0*dt)*dx);

  // caluculation of coefficient for absorption boundary condition
  // Higdon Coefficents
  a1 = 1/cos(0.0);
  a2 = 1/cos(0.0);
  d1 = 0.005;
  d2 = 0.005;
  Ca0 = (c_m0 * dt - dx)/(c_m0 * dt + dx);
  Ca1 = (dx*2) / ( c_m0 * dt + dx );
  Ca2 = (dx * c_m0 * dt * c_m0 * dt)/(dx * dx * 2 *(c_m0 * dt + dx));

  Cah1 = (a1*c_m0 * dt - dx)/(a1*c_m0 * dt + dx);
  Cah2 = (a2*c_m0 * dt - dx)/(a2*c_m0 * dt + dx);
}

// exitation condition
// pulse wave
float Driver(int T){
  float so, w, Tr, tau, al;
  w = 2 * PI * f;
  Tr = w * dt * (float)T / WN;
  switch(SC){
    case 0:
      if(Tr <= PI)so = ((1 - cos((float)(w * dt * (double)T/WN )))*sin((double)(w * dt * (double)T))/2);
      else so = sin((float)(w * dt * (double)T));
      break;

    case 1: 
      tau = (float)WN/f;
      al = (4/tau)*(4/tau);
      so = exp(-al*(dt*(float)T-tau)*(dt*(double)T-tau))*sin(w*(dt*(double)T-tau));
      break;

    case 2: 
      if(T <= W_end)so = (float)((1-cos((double)(w*dt*(double)T/WN)))*sin((double)(w * dt * (double)T))/2);
      else so = 0.0;
      break;

    case 3: //
      break;

    default: 
      tau = (float)WN/f;
      al = (4/tau)*(4/tau);
      so = exp(-al*(dt*(float)T-tau)*(dt*(double)T-tau))*sin(w*(dt*(double)T-tau));
  }
  return  so;
}

void WaveT(int id){
  if(id==3)amp1[T]= P1[SX][SY];
  else if(id==6)amp2[T]= P1[RX][RY];
}

void OutputF(void){
  char file11[50] = "field.dat";
  FILE *fp11;
  float bun;

  if(NULL == (fp11 = fopen(file11,"w"))){
    printf("\n\n Cannot Open File : %s\n",file11);
    exit(1);
  }

  for(j = 0; j <= y_max-1; j++){
    for(i = 0; i <= x_max-1; i++){
      bun = P1[i][j]*P1[i][j];
      if( bun < 1e-15)
        bun = -150;
      else 
        bun =  10.0*log10( bun );
      fprintf(fp11,"%d %d %f \t\n",i,j,(float)bun);
      i=i+mg;
    }
    fprintf(fp11,"\n");
    j=j+ng;
  }
  fclose(fp11);
  printf("output sound pressure field file : %s\n",file11);
}

void OutputW(int id){
  char file12[50] = "t_wave.dat";
  char file13[50] = "r_wave.dat";
  FILE *fp12,*fp13;
  float sub_t;
  float bun1,bun2;

  if(id==3){
    if(NULL == (fp12 = fopen(file12,"w"))){
      printf("\n\n Cannot Open File : %s\n",file12);
      exit(1);
    }
  } else if(id==6){
    if(NULL == (fp13 = fopen(file13,"w"))){
      printf("\n\n Cannot Open File : %s\n",file13);
      exit(1);
    }
  }

  for(i = 0; i <= T_max; i++){
    sub_t = dt*(float)i;
    if(id==3){
      bun1 = amp1[i];
      fprintf(fp12,"%f %f\n",sub_t,bun1);
    } else if(id==6){
      bun2 = amp2[i];
      fprintf(fp13,"%f %f\n",sub_t,bun2);
    }
  }
  if(id==3){
    fclose(fp12);
    printf("output projected wave file : %s\n",file12);
  } else if(id==6){
    fclose(fp13);
    printf("output received wave file : %s\n",file13);
  }
}

void FieldSnap(float recv[][3200]){
  FILE *fp[5000];
  char fname[8];
  int k,l,m,n;
  static  char number[11]="0123456789";
  float bun;

  k = T/del_T;
  l = (int)k/100;
  m = (int)(k%100-k%10)/10;
  n = k%10;

  fname[0] = number[l];
  fname[1] = number[m];
  fname[2] = number[n];
  fname[3] = '.';
  fname[4] = 'd';
  fname[5] = 'a';
  fname[6] = 't';
  fname[7] = 0;

  if(NULL == (fp[k] = fopen(fname,"w"))){
    printf("\n\n Cannot Open File : %s\n",fname);
    exit(1);
  }

  for(j = 0; j <= y_max-1; j++){
    for(i = 0; i <= x_max -1; i++){
      bun = recv[i][j]*recv[i][j] ;
      if(bun < 1e-15)
        bun = (float)-150;
      else
        bun = (float)(10.0*log10(bun));
      fprintf(fp[k],"%d %d %f\t\n",i,j,(float)bun);
      i = i + mg;
    }
    fprintf(fp[k],"\n");
    j = j + ng ;
  }
  fclose(fp[k]);        
}


void use_popen(void){
  FILE *fp;
  char command[MAX_STRING];
  char output[MAX_STRING];
  sprintf(command, "grep VmSize /proc/%d/status", getpid());
  if((fp = popen(command, "r")) == NULL){
    return;
  }
  while(fgets(output, MAX_STRING, fp) != NULL) {
    printf("%s", output);
  }
  if(pclose(fp) == -1){
  }
}

#endif
