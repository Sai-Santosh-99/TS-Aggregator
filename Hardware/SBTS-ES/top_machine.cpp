#include "top_machine.h"
#include "hls_math.h"
#include <math.h>

//hls::stream<ap_uint<32> > hls_output_stream;
ap_uint<32> state[ui_cycles][ui_cycle_len];
ap_uint<10> p;

void seed(ap_uint<32> s) {  // init by 32 bit seed
	#pragma HLS ARRAY_PARTITION variable=state complete dim=1

  state[0][0] = s;// & 0xFFFFFFFFUL; // for > 32 bit machines
  //state[0] = s;
  for (ap_uint<10> i = 1; i < n; ++i) {
    ap_uint<10> i_subarr=i%cycles;
    ap_uint<10> i_index =i/cycles;

    ap_uint<10> i_m1_subarr=(i-1)%cycles;
    ap_uint<10> i_m1_index =(i-1)/cycles;

    state[i_subarr][i_index] = MT_f * (state[i_m1_subarr][i_m1_index] ^ (state[i_m1_subarr][i_m1_index] >> 30)) + i;
  }
}

void mtwist_core (bool init, ap_uint<32> seed_val, ap_uint<32> stream_length, ap_uint<32> *out_stream)//, uint16_t strm_len)
{
//	#pragma HLS INTERFACE axis port=out_stream depth=1

  if(init){
    seed(seed_val);
    p=0;
  }

  ap_uint<32> x;
  ap_uint<32> temp_state;
  ap_uint<32> temp_state_next;

  ap_uint<10> q=p;
  ap_uint<10> q_subarr=q%cycles;
  ap_uint<10> q_index=q/cycles;

  bool p_wrap_next= (p==n-2);
  ap_uint<10> p_next=(p+1);
  if(p==(n-1))
    p_next=0;
  ap_uint<10> p_subarr=p%cycles;
  ap_uint<10> p_index=p/cycles;

  ap_uint<10> index_cacheA;
  ap_uint<10> index_cacheA_subarr;
  ap_uint<10> index_cacheA_index;

  ap_uint<10> index_cacheC;
  ap_uint<10> index_cacheC_subarr;
  ap_uint<10> index_cacheC_index;

  ap_uint<32> state_cacheA;//=state[m];
  ap_uint<32> state_cacheB;//=state[0];
  ap_uint<32> state_cacheC;//=state[1];

  state_cacheB=state[p_subarr][p_index];

  if(p<(n-m)){
    index_cacheA=p+m;
    index_cacheC=p+1;
  } else if(p<(n-1) ){
    index_cacheA=p+m-n;
    index_cacheC=p+1;
  } else {
    index_cacheA=m-1;
    index_cacheC=0;
    //p=0;
  }
  index_cacheA_subarr=index_cacheA%cycles;
  index_cacheA_index=index_cacheA/cycles;
  index_cacheC_subarr=index_cacheC%cycles;
  index_cacheC_index=index_cacheC/cycles;

  state_cacheA=state[index_cacheA_subarr][index_cacheA_index];
  state_cacheC=state[index_cacheC_subarr][index_cacheC_index];

  temp_state_next=state_cacheA^ twiddle(state_cacheB, state_cacheC);

  ap_uint<10> index_cacheA_next=index_cacheA+1;
  if(index_cacheA==(n-1))
    index_cacheA_next=0;
  ap_uint<10> index_cacheC_next=index_cacheC+1;
  if(index_cacheC==(n-1))
    index_cacheC_next=0;

  bool index_cacheA_wrap_next=index_cacheA_next==n-1;
  bool index_cacheC_wrap_next=index_cacheC_next==n-1;


  bool stop_next=false;
  if(stream_length==0)
	stop_next=true;

  ap_uint<32> rand_index=0;

 rand_compute: while(true){//for(ap_uint<32> rand_index=0;rand_index<stream_length; rand_index++){
    #pragma HLS DEPENDENCE variable=state array inter WAR false

	#pragma HLS PIPELINE

	if(stop_next){
      break;
	}
    stop_next=(rand_index==(stream_length-1));

    rand_index++;

    temp_state=temp_state_next;
    x = state_cacheB;//state[p];

    q=p;
    p=p_next;
    p_subarr=p%cycles;
    p_index=p/cycles;

    state_cacheB=state_cacheC;

    index_cacheA=index_cacheA_next;
    index_cacheC=index_cacheC_next;

    index_cacheA_subarr=index_cacheA%cycles;
    index_cacheA_index=index_cacheA/cycles;

    index_cacheC_index=index_cacheC/cycles;

    switch(index_cacheA_subarr){
    case 0:
      state_cacheA=state[0][index_cacheA_index];
      state_cacheC=state[4][index_cacheC_index];
      state[2][q_index]=temp_state;
      break;
    case 1:
      state_cacheA=state[1][index_cacheA_index];
      state_cacheC=state[5][index_cacheC_index];
      state[3][q_index]=temp_state;
      break;
    case 2:
      state_cacheA=state[2][index_cacheA_index];
      state_cacheC=state[6][index_cacheC_index];
      state[4][q_index]=temp_state;
      break;
    case 3:
      state_cacheA=state[3][index_cacheA_index];
      state_cacheC=state[7][index_cacheC_index];
      state[5][q_index]=temp_state;
      break;
    case 4:
      state_cacheA=state[4][index_cacheA_index];
      state_cacheC=state[0][index_cacheC_index];
      state[6][q_index]=temp_state;
      break;
    case 5:
      state_cacheA=state[5][index_cacheA_index];
      state_cacheC=state[1][index_cacheC_index];
      state[7][q_index]=temp_state;
      break;
    case 6:
      state_cacheA=state[6][index_cacheA_index];
      state_cacheC=state[2][index_cacheC_index];
      state[0][q_index]=temp_state;
      break;
    case 7:
      state_cacheA=state[7][index_cacheA_index];
      state_cacheC=state[3][index_cacheC_index];
      state[1][q_index]=temp_state;
      break;
    }

    temp_state_next=state_cacheA^ twiddle(state_cacheB, state_cacheC);


    q_subarr=p_subarr;
    q_index=p_index;

    if(p_wrap_next){
      p_next=0;
    } else {
      p_next=p+1;
    }

    p_wrap_next= (p==n-2);

    if(index_cacheA_wrap_next){
      index_cacheA_next=0;
    } else {
      index_cacheA_next++;
    }
    index_cacheA_wrap_next= (index_cacheA==n-2);


    if(index_cacheC_wrap_next){
      index_cacheC_next=0;
    } else {
      index_cacheC_next++;
    }
    index_cacheC_wrap_next= (index_cacheC==n-2);

    x ^= (x >> MT_u);
    x ^= (x << MT_s) & MT_b;
    x ^= (x << MT_t) & MT_c;
    *out_stream = x ^ (x >> MT_l);
  }
}
float thompsonQ(ap_uint<32> seed_val, ap_uint<32> X, ap_uint<32> T, ap_uint<32> N){

	float Q = 0, a = 0;
	ap_uint<32> out_rand= 0;
	unsigned long ui_core_rand = 0;
	int count = 0, index = -1;
	unsigned int ind = 0;

	unsigned int bin_count[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	bool x = false;

	for(int i=0; i<T; i++){
#pragma HLS unroll factor=64
		if(N == 4 && i == 0){
			x = true;
		}
		else{
			x = false;
		}
		mtwist_core(x, seed_val, 1, &out_rand);

		ui_core_rand= (unsigned long) out_rand;
		a = (float) ui_core_rand/ULONG_MAX;

		    if(a>=0.95)
		      ind = 19;
		    else if(a>=0.9)
		      ind = 18;
		    else if(a>=0.85)
		      ind = 17;
		    else if(a>=0.8)
		      ind = 16;
		    else if(a>=0.75)
		      ind = 15;
		    else if(a>=0.7)
		      ind = 14;
		    else if(a>=0.65)
		      ind = 13;
		    else if(a>=0.6)
		      ind = 12;
		    else if(a>=0.55)
		      ind = 11;
		    else if(a>=0.5)
		      ind = 10;
		    else if(a>=0.45)
		      ind = 9;
		    else if(a>=0.4)
		      ind = 8;
		    else if(a>=0.35)
		      ind = 7;
		    else if(a>=0.3)
		      ind = 6;
		    else if(a>=0.25)
		      ind = 5;
		    else if(a>=0.2)
		      ind = 4;
		    else if(a>=0.15)
		      ind = 3;
		    else if(a>=0.1)
		      ind = 2;
		    else if(a>=0.05)
		      ind = 1;
		    else
		      ind = 0;

		  bin_count[ind] = bin_count[ind] + 1;
	}

    final:while(count < X){
	#pragma HLS PIPELINE II=1
        index = index + 1;
        count = count + bin_count[index];
    }

    if(index==0)
      Q = 0.025;
    else if(index==1)
          Q = 0.075;
    else if(index==2)
          Q = 0.125;
    else if(index==3)
          Q = 0.175;
    else if(index==4)
          Q = 0.225;
    else if(index==5)
          Q = 0.275;
    else if(index==6)
          Q = 0.325;
    else if(index==7)
          Q = 0.375;
    else if(index==8)
          Q = 0.425;
    else if(index==9)
          Q = 0.475;
    else if(index==10)
          Q = 0.525;
    else if(index==11)
          Q = 0.575;
    else if(index==12)
          Q = 0.625;
    else if(index==13)
          Q = 0.675;
    else if(index==14)
          Q = 0.725;
    else if(index==15)
          Q = 0.775;
    else if(index==16)
          Q = 0.825;
    else if(index==17)
          Q = 0.875;
    else if(index==18)
          Q = 0.925;
    else
          Q = 0.975;

    return Q;
}

void machine(ap_uint<32> machine_val, ap_uint<32> seed_val, ap_uint<32> inform, float *Q){
#pragma HLS INTERFACE s_axilite port=inform
#pragma HLS INTERFACE axis register both port=Q

	static ap_uint<32> X, T, N;

	if(inform == 0){
		X=1;
		T=1;
		N=4;
	}
	else{
		if(inform.range(2,0) == machine_val){
			X = X + inform[3];
			T = T + 1;
		}
		N = N + 1;
	}

	*Q = thompsonQ(seed_val, X, T, N);
}

