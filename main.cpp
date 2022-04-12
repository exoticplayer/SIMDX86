#include<arm_neon.h>
#include<stdio.h>
#include<iostream>
#include<time.h>
using namespace std;
int n=50;
float **B;
void m_reset(int n)
{

    B=new float*[n];
    for(int i=0;i<n;i++)
        B[i]=new float[n];
  for(int i=0;i<n;i++)
  {
      for(int j=0;j<i;j++)
        B[i][j]=0;
      B[i][i]=1.0;
      for(int j=i+1;j<n;j++)
      {
         B[i][j]=rand();
      }

  }
  for(int k=0;k<n;k++)
  {

      for(int i=k+1;i<n;i++)
      {
          for(int j=0;j<n;j++)
            B[i][j]+=B[k][j];
      }
  }

}
void chuanxing(float **A)
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    for (int k = 0; k < n; k++)
    {
		for (int j = k + 1; j < n; j++)
		{
			A[k][j] /= A[k][k];
		}
		A[k][k] = 1.0;
		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
			{
				A[i][j] -= A[i][k] * A[k][j];
			}
			A[i][k] = 0;
		}
	}
	timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   //fprintf(result,"%lld.%09llds\n",dsec,dnsec);
   printf("%lld.%09llds",dsec,dnsec);
}
void two_paerll(float **A)
{
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);
    for(int k=0;k<n;k++)
    {
        float32x4_t vt = vdupq_n_f32(A[k][k]); // 存储的四个 float32 都初始化为A[k][k]
        int j;
        for(j=k+1;j+4<=n;j+=4)
        {
            float32x4_t va=vld1q_f32(&A[k][j]);
            va=vdivq_f32(va,vt);
            vst1q_f32(&A[k][j], va);// 将 q0 中 4 个 float32，赋值给以 d1 为起始地址的 4 个 float32

        }
        while(j<n)
            {
                A[k][j]=A[k][j]/A[k][k];
                j++;
            }
        A[k][k]=1.0;
        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
                A[i][j]-=A[i][k]*A[k][j];
            A[i][k]=0;
        }

    }
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
   //cout<<dnsec;
//fprintf(result,"%lld.%09llds\n",dsec,dnsec);

}
void three_paerll(float **A)
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    for(int k=0;k<n;k++)
    {
        for(int j=k+1;j<n;j++)
            A[k][j]/=A[k][k];
        A[k][k]=1.0;
        for(int i=k+1;i<n;i++)
        {
            float32x4_t vaik=vdupq_n_f32(A[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&A[k][j]);
                float32x4_t vaij=vld1q_f32(&A[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&A[i][j], vaij);

            }
            while(j<n)
            {
                A[i][j]-=A[k][j]*A[i][k];
                j++;
            }
            A[i][k]=0;

        }

    }
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
   //cout<<dnsec;
    //fprintf(result,"%lld.%09llds\n",dsec,dnsec);

}
void all_paerll(float **A)
{
    struct timespec sts,ets;
    timespec_get(&sts,TIME_UTC);
    for(int k=0;k<n;k++)
    {
        float32x4_t vt = vdupq_n_f32(A[k][k]); // 存储的四个 float32 都初始化为A[k][k]
        int j;
        for(j=k+1;j+4<=n;j+=4)
        {
            float32x4_t va=vld1q_f32(&A[k][j]);
            va=vdivq_f32(va,vt);
            vst1q_f32(&A[k][j], va);// 将 q0 中 4 个 float32，赋值给以 d1 为起始地址的 4 个 float32

        }
        while(j<n)
            {
                A[k][j]=A[k][j]/A[k][k];
                j++;
            }
        A[k][k]=1.0;
       for(int i=k+1;i<n;i++)
        {
            float32x4_t vaik=vdupq_n_f32(A[i][k]);
            int j=k+1;
            for(j=k+1;j+4<=n;j+=4)
            {
                float32x4_t vakj=vld1q_f32(&A[k][j]);
                float32x4_t vaij=vld1q_f32(&A[i][j]);
                float32x4_t vx=vmulq_f32(vakj,vaik);
                vaij=vsubq_f32(vaij,vx);

                vst1q_f32(&A[i][j], vaij);

            }
            while(j<n)
            {
                A[i][j]-=A[k][j]*A[i][k];
                j++;
            }
            A[i][k]=0;

        }

    }
    timespec_get(&ets,TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0)
    {
        dsec--;
        dnsec+=1000000000ll;
    }
   printf("%lld.%09llds",dsec,dnsec);
   //cout<<dnsec;
   //fprintf(result,"%lld.%09llds\n",dsec,dnsec);
    //result << (tail-head)*1000.0 / freq<<"      ";

}
//void all_paerll_duiqi()//全并行的对齐
//{
//    struct timespec sts,ets;
//    timespec_get(&sts,TIME_UTC);
//    for(int k=0;k<n;k++)
//    {
//         float32x4_t vt = vdupq_n_f32(A[k][k]);  // 存储的四个 float32 都初始化为A[k][k]
//        int j=k+1;
//        while(j%4!=0)
//        {
//            A[k][j] /= A[k][k];
//            j++;
//        }
//        for(;j+4<=n;j+=4)
//        {
//            float32x4_t va=vld1q_f32(&A[k][j]);
//            va=vdivq_f32(va,vt);
//            vst1q_f32(&A[k][j], va);// 将 q0 中 4 个 float32，赋值给以 d1 为起始地址的 4 个 float32
//
//        }
//        while(j<n)
//            {
//                A[k][j]=A[k][j]/A[k][k];
//                j++;
//            }
//        A[k][k]=1.0;
//      for(int i=k+1;i<n;i++)
//        {
//            float32x4_t vaik=vdupq_n_f32(A[i][k]);
//            int j=k+1;
//            while(j%4!=0)
//            {
//                A[i][j] -= A[i][k] * A[k][j];
//                j++;
//            }
//            for(;j+4<=n;j+=4)
//            {
//                float32x4_t vakj=vld1q_f32(&A[k][j]);
//                float32x4_t vaij=vld1q_f32(&A[i][j]);
//                float32x4_t vx=vmulq_f32(vakj,vaik);
//                vaij=vsubq_f32(vaij,vx);
//
//                vst1q_f32(&A[i][j], vaij);
//
//            }
//            while(j<n)
//            {
//                A[i][j]-=A[k][j]*A[i][k];
//                j++;
//            }
//            A[i][k]=0;
//
//        }
//
//    }
//    timespec_get(&ets,TIME_UTC);
//    time_t dsec=ets.tv_sec-sts.tv_sec;
//    long dnsec=ets.tv_nsec-sts.tv_nsec;
//    if(dnsec<0)
//    {
//        dsec--;
//        dnsec+=1000000000ll;
//    }
//   printf("%lld.%09llds",dsec,dnsec);
//
//}
int main()
{

    while(n<1500)
    {
        m_reset(n);
        float **A=new float *[n];
        for(int i=0;i<n;i++)
            A[i]=new float[n];
        for(int i=0;i<n;i++)
        {

            for(int j=0;j<n;j++)
                A[i][j]=B[i][j];
        }
        cout<<"n="<<n<<":   ";
        chuanxing(A);

        cout<<"     ";
        two_paerll(A);
        cout<<"     ";
        three_paerll(A);
        cout<<"     ";
        all_paerll(A);
        cout<<"     ";
        //all_paerll_duiqi();
        cout<<endl;
        n+=100;
    }
}
