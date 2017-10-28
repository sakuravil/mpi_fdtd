#ifndef BOUNDARY_H 
#define BOUNDARY_H

// Caution
// 
// Only Mur_1st_ABC function is parallelized.
// This function is assumed to be used only.

void Hig_2nd_ABC(void){
  for(j = 1; j <= y_max - 2; ++j){
    i=0;
    P1[i][j] = (Cah1+Cah2)*(P1[i+1][j] - PX1[0][j])
          -Cah1*Cah2*(P1[i+2][j] -2*PX1[1][j] +PX2[0][j])
          -(Cah1*(1-d2)+Cah2*(1-d1))*(PX1[2][j] -PX2[1][j])
          + ((1-d1) + (1-d2))*(PX1[1][j])
          -(1-d1)*(1-d2)*(PX2[2][j]);
    i = x_max - 1;
    P1[i][j] = (Cah1+Cah2)*(P1[i-1][j] - PX1[3][j])
          -Cah1*Cah2*(P1[i-2][j] -2*PX1[4][j] +PX2[3][j])
          -(Cah1*(1-d2)+Cah2*(1-d1))*(PX1[5][j] -PX2[4][j])
          + ((1-d1) + (1-d2))*(PX1[4][j])
          -(1-d1)*(1-d2)*(PX2[5][j]);
  }
  for(i = 1; i <= x_max - 2; ++i){
    j=0;
    P1[i][j] = (Cah1+Cah2)*(P1[i][j+1] - PY1[0][i])
          -Cah1*Cah2*(P1[i][j+2] -2*PY1[1][i] +PY2[0][i])
          -(Cah1*(1-d2)+Cah2*(1-d1))*(PY1[2][i] -PY2[1][i])
          + ((1-d1) + (1-d2))*(PY1[1][i])
          -(1-d1)*(1-d2)*(PY2[2][i]);
    j = y_max - 1;
    P1[i][j] = (Cah1+Cah2)*(P1[i][j-1] - PY1[3][i])
          -Cah1*Cah2*(P1[i][j-2] -2*PY1[4][i] +PY2[3][i])
          -(Cah1*(1-d2)+Cah2*(1-d1))*(PY1[5][i] -PY2[4][i])
          + ((1-d1) + (1-d2))*(PY1[4][i])
          -(1-d1)*(1-d2)*(PY2[5][i]);
  }
}

void Hig_1st_ABC(void){
  for(j = 1; j <= y_max - 2; ++j){
    i=0;
    P1[i][j] = (1.0-d1)*PX1[1][j] + Cah1 * (P1[i+1][j] - PX1[0][j]);
    i = x_max - 1;
    P1[i][j] = (1.0-d1)*PX1[4][j] + Cah1 * (P1[i-1][j] - PX1[3][j]);
    }
    for( i = 1 ; i <= x_max - 2 ; ++i )	{
    j=0;
    P1[i][j] = (1.0-d1)*PY1[1][i] + Cah1 * (P1[i][j+1] - PY1[0][i]);
    j = y_max - 1;
    P1[i][j] = (1.0-d1)*PY1[4][i] + Cah1 * (P1[i][j-1] - PY1[3][i]);
  }
}

void Mur_1st_ABC(int id){
  start=(x_max/P)*id;
  for(j = 1; j <= y_max - 2; ++j){
    if(id==0){
      i=0;
      P1[i][j] = PX1[1][j] + Ca0 * (P1[i+1][j] - PX1[0][j]);
    } else if(id==(P-1)){
      i = x_max - 1;
      P1[i][j] = PX1[4][j] + Ca0 * (P1[i-1][j] - PX1[3][j]);
    }
  }
  for(i = 1 + start; i <= x_max/P + start; ++i){
    if(i<=x_max-2){
      j=0;
      P1[i][j] = PY1[1][i] + Ca0 * (P1[i][j+1] - PY1[0][i]);
      j = y_max - 1;
      P1[i][j] = PY1[4][i] + Ca0 * (P1[i][j-1] - PY1[3][i]);
    }
  }
}


void Mur_2nd_ABC(void){
  float	p_1, p_2, p_3;
  i=0;
  j=1;
  P1[i][j] = PX1[1][j] + Ca0 * (P1[i+1][j] - PX1[0][j]);
  j = y_max-2;
  P1[i][j] = PX1[1][j] + Ca0 * (P1[i+1][j] - PX1[0][j]);
  i = x_max - 1;
  j = 1;
  P1[i][j] = PX1[4][j] + Ca0 * (P1[i-1][j] - PX1[3][j]);
  j = y_max-2;
  P1[i][j] = PX1[4][j] + Ca0 * (P1[i-1][j] - PX1[3][j]);
  for(j = 2; j <= y_max - 3; ++j){
    i=0;
    p_1 = Ca0*(P1[i+1][j]+PX2[0][j]);
    p_2 = Ca1*(PX1[0][j]+PX1[1][j]);
    p_3 = Ca2*(PX1[0][j+1]-2*PX1[0][j]+PX1[0][j-1]+PX1[1][j+1]-2*PX1[1][j]+PX1[1][j-1]);
    P1[i][j] = - PX2[1][j]+(p_1 + p_2 + p_3);

    i = x_max - 1;
    p_1 = Ca0*(P1[i-1][j]+PX2[3][j]);
    p_2 = Ca1*(PX1[3][j]+PX1[4][j]);
    p_3 = Ca2*(PX1[3][j+1]-2*PX1[3][j]+PX1[3][j-1]+PX1[4][j+1]-2*PX1[4][j]+PX1[4][j-1]);
    P1[i][j] = - PX2[4][j]+(p_1 + p_2 + p_3);
  }
  j=0;
  i=1;
  P1[i][j] = PY1[1][i] + Ca0 * (P1[i][j+1] - PY1[0][i]);
  i=x_max-2;
  P1[i][j] = PY1[1][i] + Ca0 * (P1[i][j+1] - PY1[0][i]);
  j = y_max - 1;
  i=1;
  P1[i][j] = PY1[4][i] + Ca0 * (P1[i][j-1] - PY1[3][i]);
  i=x_max-2;
  P1[i][j] = PY1[4][i] + Ca0 * (P1[i][j-1] - PY1[3][i]);

  for(i = 2; i <= x_max - 3; ++i){
    j=0;
    p_1 = Ca0*(P1[i][j+1]+PY2[0][i]);
    p_2 = Ca1*(PY1[0][i]+PY1[1][i]);
    p_3 = Ca2*(PY1[0][i+1]-2*PY1[0][i]+PY1[0][i-1]+PY1[1][i+1]-2*PY1[1][i]+PY1[1][i-1]);
    P1[i][j] = - PY2[1][i]+(p_1 + p_2 + p_3 );

    j = y_max - 1;
    p_1 = Ca0*(P1[i][j-1]+PY2[3][i]);
    p_2 = Ca1*(PY1[3][i]+PY1[4][i]);
    p_3 = Ca2*(PY1[3][i+1]-2*PY1[3][i]+PY1[3][i-1]+PY1[4][i+1]-2*PY1[4][i]+PY1[4][i-1]);
    P1[i][j] = - PY2[4][i]+( p_1 + p_2 + p_3);
  }
}


void save_bound(int id){
  start=(x_max/P)*id;

  for(j = 1; j <= y_max - 2; ++j){
    if(id==0){
      i=0;
      PX2[0][j] = PX1[0][j];
      PX2[1][j] = PX1[1][j];
      PX2[2][j] = PX1[2][j];
      PX1[0][j] = P1[i][j];
      PX1[1][j] = P1[i+1][j];
      PX1[2][j] = P1[i+2][j];
    } else if(id==(P-1)){
      i = x_max - 1;
      PX2[3][j] = PX1[3][j];
      PX2[4][j] = PX1[4][j];
      PX2[5][j] = PX1[5][j];
      PX1[3][j] = P1[i][j];
      PX1[4][j] = P1[i-1][j];
      PX1[5][j] = P1[i-2][j];
    }
  }
  for(i = 1 + start; i <= x_max/P + start; ++i){
    if(i<=x_max-2){
      j=0;
      PY2[0][i] = PY1[0][i];
      PY2[1][i] = PY1[1][i];
      PY2[2][i] = PY1[2][i];
      PY1[0][i] = P1[i][j];
      PY1[1][i] = P1[i][j+1];
      PY1[2][i] = P1[i][j+2];

      j = y_max - 1;
      PY2[3][i] = PY1[3][i];
      PY2[4][i] = PY1[4][i];
      PY2[5][i] = PY1[5][i];
      PY1[3][i] = P1[i][j];
      PY1[4][i] = P1[i][j-1];
      PY1[5][i] = P1[i][j-2];
    }
  }
}

void boundary(int id){
  switch(BC){
    case 0:
      Mur_1st_ABC(id); // only use
      save_bound(id);
      break;

    case 1:
      Mur_2nd_ABC();
      save_bound(id);
      break;

    case 2:
      Hig_1st_ABC();
      save_bound(id);
      break;

    case 3:
      Hig_2nd_ABC();
      save_bound(id);
      break;

    default:
      Hig_2nd_ABC();
      save_bound(id);
      break;
  }
}

#endif
