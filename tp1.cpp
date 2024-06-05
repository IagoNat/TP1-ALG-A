#include <iostream>
#include <set>
#include <gmpxx.h>
#include <vector>
#include <cmath>
#include <chrono>
#include <map>
#include "primosbulk.hpp"
using namespace std;
using namespace std::chrono;


bool miller_Rabin(const mpz_class& n, const mpz_class& a) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0) return false;

    mpz_class d = n - 1;
    int s = 0;
    while (d % 2 == 0) {
        d /= 2;
        s++;
    }

    mpz_class x = mpz_class(1);
    mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

    if(x==1 || x==n-1)
        return true;

    bool isComposite = true;
    for (int r = 1; r <= s; r++) {
        mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());
        if (x == 1) return false;
        if (x == n - 1) {
            isComposite = false;
            break;
        }
    }

    if (isComposite) return false;
    return true;
}

pair<mpz_class, mpz_class> findNextPrime(const mpz_class& N, const mpz_class& a) {
    mpz_class counter{1};
    mpz_class nextPrime = N + (N % 2 == 0 ? 1 : 2); 

    counter++;
    while (!miller_Rabin(nextPrime, a)) {
        nextPrime += 2;
        counter++;
    }
    return make_pair(counter, nextPrime);
}


pair<mpz_class, map<mpz_class, int>> factorization(mpz_class p) {
    map<mpz_class, int> factors;

    for(const mpz_class& prime: not_so_small_set_of_small_primes) {
        while((p%prime)==0) {
            p /= prime;
            factors[prime]++;
        }
    }

    if(p != 1) {
        mpz_class n_p = *(--not_so_small_set_of_small_primes.end());
        auto start = steady_clock::now();
        auto end_time = start + seconds(600);
        
        while(p != 1) {
            auto now = steady_clock::now();
            n_p = findNextPrime(n_p, 2).second;
            while((p%n_p)==0) {
                p /= n_p;
                factors[n_p]++;
            }

            if(miller_Rabin(p, 2)) {
                factors[p]++;
                p = 1;
                break;
            }
            if(now >= end_time)
                break;
        }
    }

    return make_pair(p, factors);
}

mpz_class find_generator(mpz_class p, map<mpz_class, int>& f) {
    mpz_class result {0};
    mpz_class p_minus_one = p-1;
    bool is_generator;
    for(mpz_class i=2; i<p_minus_one; i++) {
        is_generator=true;
        for(const auto& [factor, n]: f) {
            mpz_class exponent = (p_minus_one/factor);
            mpz_powm(result.get_mpz_t(), i.get_mpz_t(), exponent.get_mpz_t(), p.get_mpz_t()); 
            if(result == 1) {
                is_generator=false;
                break;
            }
        }
        if(is_generator)
            return i;
    }
    return -1;
}

mpz_class find_element_with_high_order(mpz_class p, map<mpz_class, int>& f) {
    vector<mpz_class> partial_generators;
    mpz_class p_minus_one = p-1;

    for(mpz_class i=2; i<p_minus_one && partial_generators.size() < f.size(); i++) {
        bool is_new = true;
        for(const auto& [factor, n]: f) {
            mpz_class exp = p_minus_one / factor;
            mpz_class result;
            mpz_powm(result.get_mpz_t(), i.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
            if(result == 1) {
                is_new = false;
                break;
            }
        }
        if(is_new) 
            partial_generators.push_back(i);
    }

    mpz_class element = 1;
    for(const auto& b: partial_generators) {
        element = (element * b)%p;
    }

    return element;
}

mpz_class pohlig_hellman(mpz_class g, map<mpz_class, int> f, mpz_class h, mpz_class p) {
    auto start = steady_clock::now();
    auto end_time = start + seconds(600);

    mpz_class p_minus_one = p-1;
    
    vector<mpz_class> results;
    vector<mpz_class> moduli;

    for(const auto& [factor, n]: f) {
        if(steady_clock::now() >= end_time)
                break;

        
        mpz_class f_pow_n = mpz_class(1);
        mpz_pow_ui(f_pow_n.get_mpz_t(), factor.get_mpz_t(), n);
        mpz_class g_f{g}, h_f{h};
        mpz_class f_div;
        mpz_div(f_div.get_mpz_t(), p_minus_one.get_mpz_t(), factor.get_mpz_t());
        mpz_powm(g_f.get_mpz_t(), g.get_mpz_t(), f_div.get_mpz_t(), p.get_mpz_t());

        mpz_class x_f{0};
        mpz_class f_pow_i{1};

        for(int i=0; i<n; i++) {
            mpz_class h_i, a, b;
            mpz_mul(b.get_mpz_t(), f_pow_i.get_mpz_t(), factor.get_mpz_t());
            mpz_div(a.get_mpz_t(), p_minus_one.get_mpz_t(), b.get_mpz_t());

            mpz_powm(h_i.get_mpz_t(), h_f.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
            mpz_class g_i = g_f;
            mpz_class dlog = -1;
            for(mpz_class j=0; j<factor; j++) {
                if(mpz_cmp(h_i.get_mpz_t(), g_i.get_mpz_t()) == 0) {
                    dlog = j;
                    break;
                }
                g_i = (g_i*g_f)%p;
            }
            if(dlog == -1) {
                return -1; 
            }
            x_f += dlog * f_pow_i;
            f_pow_i *= factor;
            
            mpz_class inv_f_i, s, f, w;
            mpz_div(s.get_mpz_t(), p_minus_one.get_mpz_t(), f_pow_i.get_mpz_t());
            mpz_mul(f.get_mpz_t(), dlog.get_mpz_t(), s.get_mpz_t());
            w = p_minus_one - f;
            mpz_powm(inv_f_i.get_mpz_t(), g.get_mpz_t(), w.get_mpz_t(), p.get_mpz_t());
            h_f = (h_f * inv_f_i)%p;
        }

        results.push_back(x_f);
        moduli.push_back(f_pow_n);
    }

    mpz_class result = 0;
    mpz_class product = 1;
    for(const auto& mod: moduli) {
        if(steady_clock::now() >= end_time)
                break;
        product *= mod;
    }

    for(size_t i=0; i<results.size(); i++) {
        if(steady_clock::now() >= end_time)
                break;
        mpz_class ni = product / moduli[i];
        mpz_class ni_inv;
        mpz_invert(ni_inv.get_mpz_t(), ni.get_mpz_t(), moduli[i].get_mpz_t());
        result += results[i] * ni * ni_inv;
    }

    result %= product;
    duration<double> time = steady_clock::now() - start; 
    cout << time.count() << endl;
    return result;
}

int main() {
    mpz_class N, a;
    cin >> N >> a;

    // Encontrar o menor primo p > N (determinar quantas vezes o teste de Miller-Rabin foi usado)
    pair<mpz_class, mpz_class> counter_and_prime;

    if(N <= *(--not_so_small_set_of_small_primes.end())) {
        counter_and_prime.first = 0;
        counter_and_prime.second = *not_so_small_set_of_small_primes.upper_bound(N);
    } else {
        counter_and_prime = findNextPrime(N, a);
    }

    mpz_class num_tests = counter_and_prime.first;
    mpz_class prime = counter_and_prime.second; 

    cout << prime << endl << num_tests << endl;

    // Encontrar um número g gerador de Zp. Caso seja não computável se g é gerador ou não,
    // encontrar um elemento com ordem alto e dar uma estimativa da ordem mínima.
    mpz_class g = 0;
    mpz_class h = 0;

    mpz_class p;
    map<mpz_class, int> factors;
    if(miller_Rabin(prime-1, a)) {
        factors[prime-1]++;
    } else {
        pair<mpz_class, map<mpz_class, int>> r = factorization(prime-1);
        
        p = r.first;
        factors = r.second;
    }

    if(p != 1)
        h = find_element_with_high_order(prime, factors);
    else
        g = find_generator(prime, factors);

    if(g == 0) 
        cout << h << endl;
    else
        cout << g << endl;

    // Encontrar o logaritmo discreto de a módulo p na base g. Dependendo do valor de p esta
    // resposta pode ser não computável em tempo razoável. Neste caso será necessário colocar um tempo
    // de parada. Caso tenha resposta, esta deve incluir o tempo de cálculo. 
    mpz_class result = pohlig_hellman(g, factors, a, prime);

    cout << result << endl;

    return 0;
}