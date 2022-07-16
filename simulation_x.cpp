#include "simulation_x.hpp"
#include "random_x.hpp"

XNum XSimulation::Distance(XParticle *p1, XParticle *p2) {
  XNum res = 0;
  int i;
  for (i = 0; i < D; ++i) {
    XNum a = MinimumImage1(p1->coor[i]-p2->coor[i]);
    res += a*a;
  }
  return res;
}

void XSimulation::RelativeDistance(XParticle *p1, XParticle *p2, XNum *displace) {
  int i;
  for (i = 0; i < D; ++i)
    displace[i] = MinimumImage2(p1->coor[i]-p2->coor[i]);
}

void XSimulation::Initial(std::istream &in) {
  std::string parameter;
  while (in >> parameter) {
    if (parameter == "N")
      in >> N;
    else if (parameter == "P")
      in >> P;
    else if (parameter == "T") {
      in >> T;
      beta = 1/(kB*T);
    }
    else if (parameter == "L")
      in >> L;
    else if (parameter == "seed") {
      int s;
      in >> s;
      XSetRandSeed(s);
    }
    else if (parameter == "vi")
      in >> vi;
    else if (parameter == "omg0")
      in >> omg0;
    else if (parameter == "g")
      in >> g;
    else if (parameter == "s")
      in >> s;
    else if (parameter == "h")
      in >> h;
    else if (parameter == "skip")
      in >> skip;
    else if (parameter == "step")
      in >> step;
    else if (parameter == "den_num") {
      int num;
      in >> num;
      stat.Init(num);
    }
  }
  particles = new XParticle[N*P];
  int l, j;
  XNum vs[D];
  int i;
  for (i = 0; i < D; ++i)
    vs[i] = 0;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      particles[index].m = 1;
      particles[index].n = l;
      particles[index].p = j;
      for (i = 0; i < D; ++i) {
        XNum tmp = XRandGauss() * std::sqrt(1 / (particles[index].m * beta));
        vs[i] += tmp;
        particles[index].vel[i] = tmp;
      }
      for (i = 0; i < D; ++i)
        particles[index].coor[i] = L / 2 + 1.0 * (XRandFloat()-0.5);
    }
  for (i = 0; i < D; ++i)
    vs[i] /= N*P;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (i = 0; i < D; ++i)
        particles[index].vel[i] -= vs[i];
      int j;
      for (i = 0; i < D; ++i)
        for (j = 0; j < M; ++j) {
          particles[index].theta[i][j] = 1;
          particles[index].vtheta[i][j] = 1;
          particles[index].Q[i][j] = 1.0;
        }
    }
  VelocityRescale();
  omgP = std::sqrt(P)/(beta*hBar);
  ENkCache = new XNum*[N];
  for (l = 1; l <= N; ++l)
    ENkCache[l-1] = new XNum[l];
  VBCache = new XNum[N+1];
  VBCacheE = new XNum[N+1];
  VBCacheO = new XNum[N+1];
  ForceVBCache = new XNum[D*N*P];
  ForceVBCacheE = new XNum[D*N*P];
  ForceVBCacheO = new XNum[D*N*P];
  for (i = 0; i < D; ++i)
    center.coor[i] = L/2;
  t = 0;
  ForceCache = NULL;
  count = 0;
  incre = L/2/stat.num;
  ok = true;
  even = true;
}

void XSimulation::Dump(std::ostream &out) {
  int l, j;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      out << particles[index].n << " " << particles[index].p << " " << particles[index].m << std::endl;
      for (i = 0; i < D; ++i)
        out << particles[index].coor[i] << " ";
      for (i = 0; i < D; ++i)
        out << particles[index].vel[i] << " ";
      out << std::endl;
      int j;
      for (i = 0; i < D; ++i)
        for (j = 0; j < M; ++j)
          out << particles[index].theta[i][j] << " " << particles[index].vtheta[i][j] << " " << particles[index].Q[i][j] << " ";
      out << std::endl << std::endl;
    }
}

XSimulation::~XSimulation() {
  delete[] particles;
  int l;
  for (l = 1; l <= N; ++l)
    delete[] ENkCache[l-1];
  delete[] ENkCache;
  delete[] VBCache;
  delete[] VBCacheE;
  delete[] VBCacheO;
  delete[] ForceVBCache;
  delete[] ForceVBCacheE;
  delete[] ForceVBCacheO;
  if (ForceCache)
    delete[] ForceCache;
}

void XSimulation::VelocityRescale() {
  XNum sum = 0;
  int l, j;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        sum += particles[index].m*particles[index].vel[i]*particles[index].vel[i];
    }
  XNum lam = std::sqrt((N*P)*D*T/sum);
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        particles[index].vel[i] *= lam;
    }
}

XNum XSimulation::ENk(int N2, int k) {
  XNum res = 0;
  int l, j;
  for (l = N2 - k + 1; l <= N2; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int index2 = NextIndex(l, j, N2, k);
      res += 0.5*particles[index].m*omgP*omgP*Distance(&particles[index], &particles[index2]);
    }
  return res;
}

int XSimulation::NextIndex(int l, int j, int N2, int k) {
  int res = Index(l,j+1);
  if (j == P) {
    if (l == N2)
      res = Index(N2-k+1,1);
    else
      res = Index(l+1,1);
  }
  return res;
}

int XSimulation::PrevIndex(int l, int j, int N2, int k) {
  int res = Index(l,j-1);
  if (j == 1) {
    if (l == N2 - k + 1)
      res = Index(N2,P);
    else
      res = Index(l-1,P);
  }
  return res;
}

XNum XSimulation::XExp(XNum k, XNum E, XNum EE) {
  if (vi == 0)
    return std::exp(-beta*E+EE);
  else
    return std::exp((k-1)*std::log(vi)-beta*E+EE);
}

XNum XSimulation::XMinE(int N2) {
  if (vi == 0)
    return beta*(ENkCache[N2-1][0]+VBCache[N2-1]);
  int k;
  XNum res = 100000000;
  for (k = 1; k <= N2; ++k) {
    XNum tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCache[N2-k]);
    if (tmp < res)
      res = tmp;
  }
  return res;
}

XNum XSimulation::XMinEE(int N2) {
  if (vi == 0)
    return beta*(ENkCache[N2-1][0]+VBCacheE[N2-1]);
  int k;
  XNum res = 100000000;
  for (k = 1; k <= N2; ++k) {
    XNum tmp;
    if (k%2 == 1)
      tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCacheE[N2-k]);
    else
      tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCacheO[N2-k]);
    if (tmp < res)
      res = tmp;
  }
  return res;
}

XNum XSimulation::XMinEO(int N2) {
  if (vi == 0)
    return beta*(ENkCache[N2-1][0]+VBCacheO[N2-1]);
  int k;
  XNum res = 100000000;
  for (k = 1; k <= N2; ++k) {
    XNum tmp;
    if (k%2 == 1)
      tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCacheO[N2-k]);
    else
      tmp = -(k-1)*std::log(vi)+beta*(ENkCache[N2-1][k-1]+VBCacheE[N2-k]);
    if (tmp < res)
      res = tmp;
  }
  return res;
}

void XSimulation::FillVB() {
  int N2, k;
  VBCacheE[0] = 0;
  VBCacheE[1] = ENkCache[0][0];
  VBCacheO[0] = 10000000;
  VBCacheO[1] = 10000000;
  for (N2 = 2; N2 <= N; ++N2) {
    XNum sum = 0;
    XNum tmp = XMinEE(N2);
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
      else
        sum += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
    }
    if (sum == 0)
      VBCacheE[N2] = 10000000;
    else
      VBCacheE[N2] = (tmp-std::log(sum)+std::log(N2))/beta;
    sum = 0;
    tmp = XMinEO(N2);
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
      else
        sum += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
    }
    if (sum == 0)
      VBCacheO[N2] = 10000000;
    else
      VBCacheO[N2] = (tmp-std::log(sum)+std::log(N2))/beta;
  }
  /*VBCache[0] = 0;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum sum = 0;
    XNum tmp = std::min(beta*VBCacheE[N2], beta*VBCacheO[N2]);
    sum += std::exp(-beta*VBCacheE[N2]+tmp);
    sum += std::exp(-beta*VBCacheO[N2]+tmp);
    VBCache[N2] = (tmp-std::log(sum))/beta;
  }*/
  for (N2 = 0; N2 <= N; ++N2) {
    if (even)
      VBCache[N2] = VBCacheE[N2];
    else
      VBCache[N2] = VBCacheO[N2];
  }
}

void XSimulation::dENk(int N2, int k, int l, int j, XNum *res) {
  int i;
  for (i = 0; i < D; ++i)
    res[i] = 0;
  if (l >= N2 - k + 1 && l <= N2) {
    int index = Index(l, j);
    int index2 = NextIndex(l, j, N2, k);
    int index3 = PrevIndex(l, j, N2, k);
    XNum displace[D];
    RelativeDistance(&particles[index], &particles[index2], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles[index].m*omgP*omgP*displace[i];
    RelativeDistance(&particles[index], &particles[index3], displace);
    for (i = 0; i < D; ++i)
      res[i] += particles[index].m*omgP*omgP*displace[i];
  }
}

void XSimulation::FillForceVB() {
  int N2, k, l, j;
  XNum *res = new XNum[D*(N+1)];
  XNum *res2 = new XNum[D*(N+1)];
  XNum grad[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int i;
      for (i = 0; i < D; ++i)
        res[i] = 0;
      for (i = 0; i < D; ++i)
        res2[i] = 0;
      for (N2 = 1; N2 <= N; ++N2) {
        XNum sum2 = 0;
        XNum tmp = XMinEE(N2);
        for (k = 1; k <= N2; ++k) {
          if (vi == 0 && k-1 != 0)
            continue;
          if (k%2 == 1)
            sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
          else
            sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
        }
        for (i = 0; i < D; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            if (vi == 0 && k-1 != 0)
              continue;
            dENk(N2, k, l, j, grad);
            if (k%2 == 1)
              sum += (grad[i]+res[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
            else
              sum += (grad[i]+res2[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
          }
          if (sum2 == 0)
            res[D*N2+i] = 0;
          else
            res[D*N2+i] = sum / sum2;
        }
        sum2 = 0;
        tmp = XMinEO(N2);
        for (k = 1; k <= N2; ++k) {
          if (vi == 0 && k-1 != 0)
            continue;
          if (k%2 == 1)
            sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
          else
            sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
        }
        for (i = 0; i < D; ++i) {
          XNum sum = 0;
          for (k = 1; k <= N2; ++k) {
            if (vi == 0 && k-1 != 0)
              continue;
            dENk(N2, k, l, j, grad);
            if (k%2 == 1)
              sum += (grad[i]+res2[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
            else
              sum += (grad[i]+res[D*(N2-k)+i])*XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
          }
          if (sum2 == 0)
            res2[D*N2+i] = 0;
          else
            res2[D*N2+i] = sum / sum2;
        }
      }
      int index = Index(l, j);
      for (i = 0; i < D; ++i) {
        ForceVBCacheE[D*index+i] = res[D*N+i];
        ForceVBCacheO[D*index+i] = res2[D*N+i];
      }
    }
  if (even) {
    for (int i = 0; i < N*P*D; ++i)
      ForceVBCache[i] = ForceVBCacheE[i];
  }
  else {
    for (int i = 0; i < N*P*D; ++i)
      ForceVBCache[i] = ForceVBCacheO[i];
  }
  delete[] res;
  delete[] res2;
}

void XSimulation::FillENk() {
  int l, j;
  for (l = 1; l <= N; ++l) {
    for (j = 1; j <= l; ++j)
      ENkCache[l-1][j-1] = ENk(l, j);
  }
}

void XSimulation::TrapForce(int index, XNum *res) {
  int i;
  XNum displace[D];
  RelativeDistance(&particles[index], &center, displace);
  for (i = 0; i < D; ++i)
    res[i] = -particles[index].m*omg0*omg0*displace[i];
}

void XSimulation::PairForce(int index, int index2, XNum *res) {
  int i;
  XNum displace[D];
  //XNum inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles[index], &particles[index2])/(s*s));
  XNum inter = g/std::pow(Distance(&particles[index], &particles[index2]), 1.5);
  RelativeDistance(&particles[index], &particles[index2], displace);
  for (i = 0; i < D; ++i)
    res[i] = displace[i]*inter;
    //res[i] = ((2*displace[i])/(s*s))*inter;
}

XNum *XSimulation::Force() {
  XNum *res = new XNum[D*N*P];
  FillENk();
  FillVB();
  FillForceVB();
  int l, j, k;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        res[D*index+i] = -ForceVBCache[D*index+i]/particles[index].m;
    }
  XNum trap[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      TrapForce(index, trap);
      int i;
      for (i = 0; i < D; ++i)
        res[D*index+i] += trap[i]/particles[index].m/P;
    }
  XNum inter[D];
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = Index(k, j);
        PairForce(index, index2, inter);
        int i;
        for (i = 0; i < D; ++i)
          res[D*index+i] += inter[i]/particles[index].m/P;
      }
    }
  return res;
}

void XSimulation::DumpForce(std::ostream &out) {
  int l, j;
  for (l = 1; l <= N; ++l) {
    for (j = 1; j <= l; ++j)
      out << ENkCache[l-1][j-1] << " ";
    out << std::endl;
  }
  out << std::endl;
  for (l = 0; l <= N; ++l)
    out << VBCache[l] << " ";
  out << std::endl;
  out << std::endl;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        out << ForceVBCache[D*index+i] << " ";
    }
  out << std::endl;
  out << std::endl;
  XNum *f = Force();
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      int i;
      for (i = 0; i < D; ++i)
        out << f[D*index+i] << " ";
    }
  out << std::endl;
  out << std::endl;
  delete[] f;
}

XNum XSimulation::VBEnergy() {
  XNum *res = new XNum[N+1];
  XNum *res2 = new XNum[N+1];
  res[0] = 0;
  res2[0] = 0;
  int N2, k;
  for (N2 = 1; N2 <= N; ++N2) {
    XNum tmp = XMinEE(N2);
    XNum sum2 = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
      else
        sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
    }
    XNum sum = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum += (res[N2-k]-ENkCache[N2-1][k-1])*XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
      else
        sum += (res2[N2-k]-ENkCache[N2-1][k-1])*XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
    }
    if (sum2 == 0)
      res[N2] = 0;
    else
      res[N2] = sum / sum2;
    tmp = XMinEO(N2);
    sum2 = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
      else
        sum2 += XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
    }
    sum = 0;
    for (k = 1; k <= N2; ++k) {
      if (vi == 0 && k-1 != 0)
        continue;
      if (k%2 == 1)
        sum += (res2[N2-k]-ENkCache[N2-1][k-1])*XExp(k,ENkCache[N2-1][k-1]+VBCacheO[N2-k],tmp);
      else
        sum += (res[N2-k]-ENkCache[N2-1][k-1])*XExp(k,ENkCache[N2-1][k-1]+VBCacheE[N2-k],tmp);
    }
    if (sum2 == 0)
      res2[N2] = 0;
    else
      res2[N2] = sum / sum2;
  }
  XNum e;
  if (even)
    e = res[N];
  else
    e = res2[N];
  delete[] res;
  delete[] res2;
  return e;
}

XNum XSimulation::TrapEnergy(int index) {
  return 0.5*particles[index].m*omg0*omg0*Distance(&particles[index], &center);
}

XNum XSimulation::PairEnergy(int index, int index2) {
  //XNum inter = (g/(M_PI*s*s))*std::exp(-Distance(&particles[index], &particles[index2])/(s*s));
  XNum inter = g/std::sqrt(Distance(&particles[index], &particles[index2]));
  return 0.5*inter;
}

XNum XSimulation::TotalEnergy() {
  XNum e = 0;
  e += P*D*N/(2*beta);
  e += VBEnergy();
  int l, j, k;
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      e += TrapEnergy(index)/P;
    }
  for (l = 1; l <= N; ++l)
    for (j = 1; j <= P; ++j) {
      int index = Index(l, j);
      for (k = 1; k <= N; ++k) {
        if (l == k) continue;
        int index2 = Index(k, j);
        e += PairEnergy(index, index2)/P;
      }
    }
  return e;
}

void XSimulation::DumpEnergy(std::ostream &out) {
  out << VBEnergy() << " " << TotalEnergy() << std::endl;
  out << std::endl;
}

void XSimulation::NHForce(XNum m, int f, XNum *v, XNum *Q, XNum *vtheta, XNum *res) {
  XNum sum = 0;
  int i;
  for (i = 0; i < f; ++i)
    sum += m*v[i]*v[i];
  res[0] = (sum-f*(1/beta))/Q[0];
  for (i = 1; i < M; ++i)
    res[i] = (Q[i-1]*vtheta[i-1]*vtheta[i-1]-(1/beta))/Q[i];
}

void XSimulation::UpdateMNHC_VV3() {
  if (ForceCache == NULL)
    ForceCache = Force();
  int j;
  XNum NHF[M];
  XNum v[1];
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (IsNan(ForceCache[D*j+i]) || IsInf(ForceCache[D*j+i])) {
        ok = false;
        return;
      }
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      particles[j].vel[i] = particles[j].vel[i]*std::exp(-0.5*h*particles[j].vtheta[i][0])+0.5*h*ForceCache[D*j+i]*std::exp(-0.25*h*particles[j].vtheta[i][0]);
      int M2 = M/2;
      int k;
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-2] = particles[j].theta[i][2*k-2]+h*particles[j].vtheta[i][2*k-2]/2;
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-1] = particles[j].vtheta[i][2*k-1]*std::exp(-0.5*h*((k==M2)?0:particles[j].vtheta[i][2*k]))+0.5*h*NHF[2*k-1]*std::exp(-0.25*h*((k==M2)?0:particles[j].vtheta[i][2*k]));
      particles[j].coor[i] = particles[j].coor[i]+h*particles[j].vel[i];
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-1] = particles[j].theta[i][2*k-1]+h*particles[j].vtheta[i][2*k-1];
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-2] = particles[j].vtheta[i][2*k-2]*std::exp(-h*particles[j].vtheta[i][2*k-1])+h*NHF[2*k-2]*std::exp(-0.5*h*particles[j].vtheta[i][2*k-1]);
    }
  }
  t += h;
  delete[] ForceCache;
  ForceCache = Force();
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (IsNan(ForceCache[D*j+i]) || IsInf(ForceCache[D*j+i])) {
        ok = false;
        return;
      }
      int M2 = M/2;
      particles[j].vel[i] = particles[j].vel[i]*std::exp(-0.5*h*particles[j].vtheta[i][0])+0.5*h*ForceCache[D*j+i]*std::exp(-0.25*h*particles[j].vtheta[i][0]);
      int k;
      for (k = 1; k <= M2; ++k)
        particles[j].theta[i][2*k-2] = particles[j].theta[i][2*k-2]+h*particles[j].vtheta[i][2*k-2]/2;
      v[0] = particles[j].vel[i];
      NHForce(particles[j].m, 1, v, particles[j].Q[i], particles[j].vtheta[i], NHF);
      for (k = 1; k <= M2; ++k)
        particles[j].vtheta[i][2*k-1] = particles[j].vtheta[i][2*k-1]*std::exp(-0.5*h*((k==M2)?0:particles[j].vtheta[i][2*k]))+0.5*h*NHF[2*k-1]*std::exp(-0.25*h*((k==M2)?0:particles[j].vtheta[i][2*k]));
    }
  }
}

void XSimulation::PeriodBoundary() {
  int j;
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i) {
      if (particles[j].coor[i] < 0) {
        int n = (int)(std::abs(particles[j].coor[i]) / L);
        particles[j].coor[i] += (n + 1) * L;
      }
      else {
        int n = (int)(std::abs(particles[j].coor[i]) / L);
        particles[j].coor[i] -= n * L;
      }
    }
  }
}

XNum XSimulation::Temperature() {
  XNum sum = 0;
  int j;
  for (j = 0; j < N*P; ++j) {
    int i;
    for (i = 0; i < D; ++i)
      sum += particles[j].m*particles[j].vel[i]*particles[j].vel[i];
  }
  return sum/(D*N*P);
}

void XSimulation::MakeStat() {
  ++count;
  stat.e += TotalEnergy();
  stat.T += Temperature();
  int j;
  XNum c = 0;
  for (j = 0; j < N*P; ++j)
    c += particles[j].coor[0];
  c /= N*P;
  XNum VBN;
  if (even)
    VBN = VBCacheE[N];
  else
    VBN = VBCacheO[N];
  XNum tmpv = vi;
  vi = vi2;
  FillVB();
  for (j = 0; j < N*P; ++j) {
    //XNum dis = std::sqrt(Distance(&particles[j], &center));
    //XNum dis = std::abs(MinimumImage2(particles[j].coor[0]-c));
    XNum dis;
    incre = 100.0/stat.num;
    if (even)
      dis = -beta*(VBCacheO[N]-VBN)+50;
    else
      dis = -beta*(VBCacheE[N]-VBN)+50;
    if (dis < 0) continue;
    int index = int(dis/incre);
    if (index < 0)
      index = 0;
    if (index >= stat.num)
      continue;
    stat.den[index] += 1;
  }
  /*if (even)
    stat.sign += std::exp(-beta*(VBCacheO[N]-VBCacheE[N]));
  else
    stat.sign += std::exp(-beta*(VBCacheE[N]-VBCacheO[N]));*/
  XNum tmp;
  if (even) {
    tmp = std::min(beta*VBN, beta*VBCacheO[N]);
    stat.sign += std::exp(-beta*VBCacheO[N]+tmp)/std::exp(-beta*VBN+tmp);
  }
  else {
    tmp = std::min(beta*VBCacheE[N], beta*VBN);
    stat.sign += std::exp(-beta*VBCacheE[N]+tmp)/std::exp(-beta*VBN+tmp);
  }
  vi = tmpv;
}

void XSimulation::DumpStat(std::ostream &out) {
  out << "E " << stat.e << std::endl;
  out << "T " << stat.T << std::endl;
  out << "s " << stat.sign << std::endl;
  int i;
  for (i = 0; i < stat.num; ++i)
    out << "{" << i*incre-50 << "," << stat.den[i] << "},";
  out << std::endl;
  out << std::endl;
}

void XSimulation::NormalizeStat() {
  stat.e /= count;
  stat.T /= count;
  stat.sign /= count;
  int i;
  //for (i = 0; i < stat.num; ++i)
    //stat.den[i] /= std::pow((i+0.5)*incre,D-1);
  XNum norm = 0;
  for (i = 0; i < stat.num; ++i)
    norm += stat.den[i]*incre;
  for (i = 0; i < stat.num; ++i)
    stat.den[i] /= norm;
}