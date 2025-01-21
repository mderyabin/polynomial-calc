#include <iostream>
#include <random>
#include <cmath>

using namespace std;

void print_vector(uint64_t *ax, size_t n);
void generate_vector(uint64_t *ax, size_t n, uint64_t m);
void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m);
uint64_t powmod(uint64_t a, uint64_t e, uint64_t m);

uint64_t mod_red_F(uint64_t a, size_t K);
uint64_t leftshift_cyclic(uint64_t a, size_t shift, size_t K);

int main(int argc, char const *argv[]) {

    uint64_t K = 8;
    uint64_t F = (1ULL << K) + 1;
    uint64_t F2 = (1ULL << (2*K));
    uint64_t N = 4;
    uint64_t N_inv = (1<<(2*K - 2))%F;

    size_t w_pow = K / N;
    uint64_t w = 1ULL << w_pow;
    uint64_t w_inv = (1<<(2*K - w_pow)) % F;


    cout << "params" << endl;
    cout << "K = " << K << endl;
    cout << "F = " << F << endl;
    cout << "N = " << N << endl;
    cout << "N_inv = " << N_inv << endl;
    cout << "check = " << (N*N_inv)%F << endl;

    cout << "w = " << w << endl;
    cout << "w_inv = " << w_inv << endl;
    cout << "check = " << (w*w_inv)%F << endl;


    size_t   ss = 19;
    uint64_t  q = 56;
    uint64_t q1 = q * (1<<ss);
    uint64_t q2 = q1 % F;
    uint64_t q3 = leftshift_cyclic(q, ss, K);
    uint64_t q4 = mod_red_F(q3, K);


    cout << "q  = " << q << endl;
    cout << "q1 = " << q1 << endl;
    cout << "q2 = " << q2 << endl;
    cout << "q3 = " << q3 << endl;
    cout << "q4 = " << q4 << endl;

    uint64_t ax[N];
    uint64_t ax_fft[N];
    uint64_t bx[N];
    uint64_t bx_fft[N];
    uint64_t cx[N];
    uint64_t cx_fft[N];
    uint64_t cx_ref[N];

    generate_vector(ax, N, F);
    generate_vector(bx, N, F);

    NaiveNegacyclicConvolution(cx_ref, ax, bx, N, F);

    cout << "    ax = "; print_vector(ax, N);     cout << endl;
    cout << "    bx = "; print_vector(bx, N);     cout << endl;
    cout << "cx_ref = "; print_vector(cx_ref, N); cout << endl;

    
    for (size_t i = 0; i < N; i++) {
        ax_fft[i] = 0;
        bx_fft[i] = 0;
        for (size_t j = 0; j < N; j++) {
            // uint64_t tt = powmod(w, 2*i*j + j, F);
            // ax_fft[i] += (tt*ax[j])%F; 
            uint64_t rr = leftshift_cyclic(ax[j], (w_pow*(2*i*j + j)) % (2*K), K);
            ax_fft[i] += mod_red_F(rr, K);

            uint64_t tt = leftshift_cyclic(bx[j], (w_pow*(2*i*j + j)) % (2*K), K);
            bx_fft[i] += mod_red_F(tt, K);
        }
        ax_fft[i] %= F;
        bx_fft[i] %= F;
        cx_fft[i] = (ax_fft[i] * bx_fft[i]) % F;
    }

    cout << "ax_fft = "; print_vector(ax_fft, N); cout << endl;
    cout << "bx_fft = "; print_vector(bx_fft, N); cout << endl;

    
    for (size_t i = 0; i < N; i++) {
        cx[i] = 0;
        for (size_t j = 0; j < N; j++) {
            // uint64_t tt = powmod(w_inv, 2*i*j + i, F);
            // cx[i] += (tt*ax_fft[j])%F; 
            uint64_t rr = leftshift_cyclic(cx_fft[j], ((2*K - w_pow)*(2*i*j + i)) % (2*K), K);
            cx[i] += mod_red_F(rr, K);
        }
        cx[i] %= F;

        cx[i] = (cx[i]*N_inv) % F;
    }

    cout << "    cx = "; print_vector(cx, N); cout << endl;

    return 0;
}

void generate_vector(uint64_t *ax, size_t n, uint64_t m) {
    std::random_device rd;
    std::default_random_engine en(rd());
    std::uniform_int_distribution<int> rnd(0, m-1);

    for (size_t i = 0; i < n; i++)
        ax[i] = rnd(en);
}


void NaiveNegacyclicConvolution(uint64_t *cx, const uint64_t *ax, const uint64_t *bx, size_t N, uint64_t m) {
    for (size_t i = 0; i < N; i++) {
        cx[i] = 0;
        for (size_t j = 0; j <= i; j++) {
            cx[i] += (ax[j] * bx[i - j]) % m;
        }
        for (size_t j = i+1; j <= N-1; j++) {
            cx[i] -= (ax[j] * bx[N + i - j]) % m;
        }
        cx[i] %= m;
    }   
}

void print_vector(uint64_t *ax, size_t n) {
    cout << "[" << ax[0];
    for (size_t i = 1; i < n; i++) {
        cout << ", " << ax[i];
    }
    cout << "]";
}

uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t t;
    for(t = 1; e; e >>= 1) {
        if (e & 1)
            t = (t * a) % m;
        a = (a * a) % m;
    }
    return t;
}

uint64_t mod_red_F(uint64_t a, size_t K) {
    uint64_t mask = (1ULL << K) - 1;
    uint64_t F = (1ULL << K) + 1;
    uint64_t a_lo = a & mask;
    uint64_t a_hi = a >> K;

    uint64_t res = (a_lo > a_hi) ? a_lo - a_hi : F - a_hi + a_lo;
    return res;
}

uint64_t leftshift_cyclic(uint64_t a, size_t shift, size_t K) {
    uint64_t mask = (1ULL << (2*K))  - 1;
    a = a << shift;
    uint64_t a_hi = a >> (2*K);

    uint64_t res = (a&mask) + a_hi;

    return res;
}