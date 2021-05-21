#include <hls_math.h>
#include <math.h>
#include <ap_int.h>
#include <ap_fixed.h>

typedef ap_ufixed<27, 15> my_type;

my_type checkUCB(ap_uint<32> X, ap_uint<32> T, ap_uint<32> N){
#pragma HLS INTERFACE axis register both port=X
#pragma HLS INTERFACE axis register both port=T
#pragma HLS INTERFACE axis register both port=N
#pragma HLS INTERFACE axis register both port=Q

	my_type Q_in, a, b, c, d;

	a = (my_type) X/T;
	b = (my_type) hls::log(N);
	c = (my_type) b/T;
	d = (my_type) hls::sqrt(c);

	Q_in = (my_type) a + d;

	return Q_in;
}

void machine(ap_uint<32> machine_val, ap_uint<32> inform, my_type *Q){
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

	*Q = checkUCB(X, T, N);
}
