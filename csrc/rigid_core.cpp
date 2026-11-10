// C++ core of pydicomrt's rigid_registration(algorithm="stochastic"), loaded by
// pydicomrt/reg/method/_rigid_core.py.
//
// The algorithm is pydicomrt/reg/method/_stochastic_rigid.py; this file is that module's inner loops
// with every per-sample step fused, parallelised over a small thread pool. The two must stay in
// step -- the parity tests in test/reg/test_rigid_stochastic.py compare them sample for sample.
//
// Determinism: samples are split into a fixed number of chunks, each with its own partial histogram
// and gradient, summed in chunk order, so the result does not depend on the thread count. Random
// coordinates come from a counter-based generator: sample k of draw d is a pure function of
// (seed, d, k), which the numpy implementation reproduces bit for bit.
//
// Plain C ABI, loaded through ctypes, so one wheel per platform serves every Python version.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#ifdef _WIN32
#define EXPORT extern "C" __declspec(dllexport)
#else
#define EXPORT extern "C" __attribute__((visibility("default")))
#endif

extern "C" {

struct Vol {
    const float* v;       // (z, y, x) C order, as sitk.GetArrayFromImage
    int64_t size[3];      // x, y, z
    double origin[3];
    double to_phys[9];    // row-major: direction * spacing
    double to_index[9];   // its inverse
    double vmin, vmax;
};

struct Level {
    const Vol* f;
    const Vol* m;
    double centre[3];
    int bins;
    uint64_t seed;
    int threads;
};

}  // extern "C"

namespace {

constexpr int NCHUNK = 64;
constexpr double ASGD_A = 20.0, ASGD_ALPHA = 0.602, SIGMOID_MAX = 1.0, SIGMOID_MIN = -0.8;
constexpr int64_t RHO_SAMPLES = 20000;
// No step moves points by more than this many step_mm. The gain is set from the gradient at the
// start of a level, which is ~0 when that start is already the optimum (two identical images);
// the next, ordinary gradient times that gain would then throw the transform tens of mm away.
constexpr double MAX_STEP_RATIO = 2.0;

inline uint64_t splitmix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

inline double uniform(uint64_t seed, uint64_t draw, uint64_t k) {
    uint64_t h = splitmix(seed ^ splitmix(draw * 0xD1B54A32D192ED03ull ^ splitmix(k)));
    return (h >> 11) * (1.0 / 9007199254740992.0);
}

inline double bspline3(double u) {
    double a = std::fabs(u);
    if (a < 1) return (4 - 6 * a * a + 3 * a * a * a) / 6;
    if (a < 2) return (2 - a) * (2 - a) * (2 - a) / 6;
    return 0;
}

inline double bspline3_d(double u) {
    double a = std::fabs(u);
    if (a < 1) return -2 * u + 1.5 * u * a;
    if (a < 2) return u > 0 ? -(2 - a) * (2 - a) / 2 : (2 - a) * (2 - a) / 2;
    return 0;
}

inline void matvec(const double* M, const double* x, double* y) {
    for (int i = 0; i < 3; ++i) y[i] = M[3 * i] * x[0] + M[3 * i + 1] * x[1] + M[3 * i + 2] * x[2];
}

// Trilinear value at continuous index u (assumed inside); optionally the exact index-space gradient.
inline double trilinear(const Vol& V, const double* u, double* g) {
    int64_t i[3];
    double f[3];
    for (int d = 0; d < 3; ++d) {
        int64_t k = (int64_t)std::floor(u[d]);
        if (k < 0) k = 0;
        if (k > V.size[d] - 2) k = V.size[d] - 2;
        i[d] = k;
        f[d] = u[d] - k;
    }
    const int64_t nx = V.size[0], nxy = V.size[0] * V.size[1];
    const float* p = V.v + (i[2] * V.size[1] + i[1]) * nx + i[0];
    double c000 = p[0], c100 = p[1], c010 = p[nx], c110 = p[nx + 1];
    double c001 = p[nxy], c101 = p[nxy + 1], c011 = p[nxy + nx], c111 = p[nxy + nx + 1];
    double fx = f[0], fy = f[1], fz = f[2];
    double c00 = c000 + (c100 - c000) * fx, c10 = c010 + (c110 - c010) * fx;
    double c01 = c001 + (c101 - c001) * fx, c11 = c011 + (c111 - c011) * fx;
    double c0 = c00 + (c10 - c00) * fy, c1 = c01 + (c11 - c01) * fy;
    if (g) {
        g[0] = ((c100 - c000) * (1 - fy) + (c110 - c010) * fy) * (1 - fz) +
               ((c101 - c001) * (1 - fy) + (c111 - c011) * fy) * fz;
        g[1] = (c10 - c00) * (1 - fz) + (c11 - c01) * fz;
        g[2] = c1 - c0;
    }
    return c0 + (c1 - c0) * fz;
}

// Runs body(c) for every chunk c in [0, NCHUNK) on `threads` threads, the caller included.
//
// Not OpenMP: libgomp's default wait spins for milliseconds, and on a machine whose cores other jobs
// were using that made a 24-thread run ten times slower than an 8-thread one. Idle workers here spin
// for SPIN_US -- long enough to bridge the serial gaps between one parallel loop and the next -- then
// block. Chunks are claimed from a ticket that carries the run's generation in its upper half, so a
// worker that wakes late cannot claim a chunk of a run it did not see published.
class Pool {
public:
    explicit Pool(int threads) {
        for (int i = 1; i < threads; ++i) workers.emplace_back([this] { loop(); });
    }
    ~Pool() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            stop = true;
        }
        wake.notify_all();
        for (auto& w : workers) w.join();
    }
    void run(const std::function<void(int)>& f) {
        body = &f;
        pending.store(NCHUNK, std::memory_order_relaxed);
        uint64_t gen = (ticket.load() >> 32) + 1;
        ticket.store(gen << 32);  // publishes body and pending
        if (sleeping.load() > 0) {
            std::lock_guard<std::mutex> lock(mutex);
            wake.notify_all();
        }
        work(gen);
        if (!spin([&] { return pending.load() == 0; })) {
            std::unique_lock<std::mutex> lock(mutex);
            finished.wait(lock, [&] { return pending.load() == 0; });
        }
    }

private:
    static constexpr int SPIN_US = 200;

    template <class Ready>
    static bool spin(Ready ready) {
        auto until = std::chrono::steady_clock::now() + std::chrono::microseconds(SPIN_US);
        for (int i = 0;; ++i) {
            if (ready()) return true;
            if ((i & 63) == 63 && std::chrono::steady_clock::now() > until) return false;
        }
    }
    void work(uint64_t gen) {
        uint64_t t = ticket.load();
        while ((t >> 32) == gen && (t & 0xffffffffu) < NCHUNK) {
            if (!ticket.compare_exchange_weak(t, t + 1)) continue;
            (*body)((int)(t & 0xffffffffu));
            if (pending.fetch_sub(1) == 1) {
                std::lock_guard<std::mutex> lock(mutex);
                finished.notify_all();
            }
            t = ticket.load();
        }
    }
    void loop() {
        uint64_t seen = 0;
        for (;;) {
            auto published = [&] { return (ticket.load() >> 32) != seen; };
            if (!spin(published)) {
                std::unique_lock<std::mutex> lock(mutex);
                ++sleeping;
                wake.wait(lock, [&] { return stop || published(); });
                --sleeping;
                if (stop) return;
            }
            seen = ticket.load() >> 32;
            work(seen);
        }
    }
    std::vector<std::thread> workers;
    std::mutex mutex;
    std::condition_variable wake, finished;
    const std::function<void(int)>* body = nullptr;
    std::atomic<uint64_t> ticket{0};
    std::atomic<int> pending{0}, sleeping{0};
    bool stop = false;
};

struct Sample {
    int fb;
    double mb;
    double dw[3], dt[3];
};

class Problem {
public:
    explicit Problem(const Level& L) : L(L), F(*L.f), M(*L.m), B(L.bins), pool(std::max(1, L.threads)) {
        fscale = (B - 1) / std::max(F.vmax - F.vmin, 1e-12);
        mscale = (B - 5) / std::max(M.vmax - M.vmin, 1e-12);
        hist.assign((size_t)NCHUNK * B * B, 0.0);
        part.assign((size_t)NCHUNK * 6, 0.0);
        count.assign(NCHUNK, 0);
    }

    // Fixed sample k of draw d: continuous index u and physical point minus the centre.
    void fixed_point(uint64_t d, uint64_t k, double* u, double* xc) const {
        for (int a = 0; a < 3; ++a) u[a] = uniform(L.seed, d, 3 * k + a) * (F.size[a] - 1);
        matvec(F.to_phys, u, xc);
        for (int a = 0; a < 3; ++a) xc[a] += F.origin[a] - L.centre[a];
    }

    // RMS distance of the first RHO_SAMPLES points of draw d from the centre: the rotation scale.
    double rho(uint64_t d) const {
        double s = 0;
        for (int64_t k = 0; k < RHO_SAMPLES; ++k) {
            double u[3], xc[3];
            fixed_point(d, k, u, xc);
            s += xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2];
        }
        return std::sqrt(s / RHO_SAMPLES);
    }

    // MI and its gradient (dw, dt) on n samples of draw d, or on the given fixed indices u_given.
    // Returns the number of samples that map inside the moving image; below 64 nothing is written.
    int64_t evaluate(const double* R, const double* t, uint64_t d, int64_t n, const double* u_given,
                     double* mi_out, double* grad_out) {
        buf.resize(n);
        inside.resize(n);
        std::fill(hist.begin(), hist.end(), 0.0);
        // Per-chunk sums stay in locals until the chunk ends: neighbouring chunks' slots share
        // cache lines, and chunks are not handed to threads in contiguous runs.
        pool.run([&](int c) {
            double* h = hist.data() + (size_t)c * B * B;
            int64_t lo = n * c / NCHUNK, hi = n * (c + 1) / NCHUNK, cnt = 0;
            for (int64_t k = lo; k < hi; ++k) {
                double u[3], xc[3], r[3], y[3], um[3];
                if (u_given) {
                    for (int a = 0; a < 3; ++a) u[a] = u_given[3 * k + a];
                    matvec(F.to_phys, u, xc);
                    for (int a = 0; a < 3; ++a) xc[a] += F.origin[a] - L.centre[a];
                } else {
                    fixed_point(d, k, u, xc);
                }
                matvec(R, xc, r);
                for (int a = 0; a < 3; ++a) y[a] = r[a] + L.centre[a] + t[a] - M.origin[a];
                matvec(M.to_index, y, um);
                bool in = true;
                for (int a = 0; a < 3; ++a) in = in && um[a] >= 0 && um[a] <= M.size[a] - 1;
                inside[k] = in;
                if (!in) continue;
                double fv = trilinear(F, u, nullptr);
                double gi[3];
                double mv = trilinear(M, um, gi);
                int fb = (int)std::floor((fv - F.vmin) * fscale + 0.5);
                fb = fb < 0 ? 0 : (fb > B - 1 ? B - 1 : fb);
                double mb = (mv - M.vmin) * mscale + 2.0;
                mb = mb < 2.0 ? 2.0 : (mb > B - 3.0 ? B - 3.0 : mb);
                int tap0 = (int)std::floor(mb) - 1;
                for (int j = 0; j < 4; ++j) h[fb * B + tap0 + j] += bspline3(mb - (tap0 + j));
                // d moving / d y: g_phys[j] = sum_i gi[i] to_index[i][j]
                double gp[3];
                for (int j = 0; j < 3; ++j)
                    gp[j] = gi[0] * M.to_index[j] + gi[1] * M.to_index[3 + j] + gi[2] * M.to_index[6 + j];
                Sample& s = buf[k];
                s.fb = fb;
                s.mb = mb;
                s.dw[0] = (r[1] * gp[2] - r[2] * gp[1]) * mscale;
                s.dw[1] = (r[2] * gp[0] - r[0] * gp[2]) * mscale;
                s.dw[2] = (r[0] * gp[1] - r[1] * gp[0]) * mscale;
                for (int j = 0; j < 3; ++j) s.dt[j] = gp[j] * mscale;
                cnt++;
            }
            count[c] = cnt;
        });
        int64_t total = 0;
        for (int c = 0; c < NCHUNK; ++c) total += count[c];
        if (total < 64) return total;

        std::vector<double> joint((size_t)B * B, 0.0);
        for (int c = 0; c < NCHUNK; ++c) {
            const double* h = hist.data() + (size_t)c * B * B;
            for (int i = 0; i < B * B; ++i) joint[i] += h[i];
        }
        double sum = 0;
        for (double v : joint) sum += v;
        std::vector<double> pf(B, 0.0), pm(B, 0.0);
        for (int i = 0; i < B; ++i)
            for (int j = 0; j < B; ++j) {
                double p = joint[i * B + j] /= sum;
                pf[i] += p;
                pm[j] += p;
            }
        double mi = 0;
        Lmat.assign((size_t)B * B, 0.0);
        for (int i = 0; i < B; ++i)
            for (int j = 0; j < B; ++j) {
                double p = joint[i * B + j];
                if (p > 0) {
                    mi += p * (std::log(p) - std::log(pf[i] * pm[j] + 1e-300));
                    Lmat[i * B + j] = std::log(p) - std::log(pm[j] + 1e-300);
                }
            }
        *mi_out = mi;
        if (!grad_out) return total;

        pool.run([&](int c) {
            double g[6] = {0, 0, 0, 0, 0, 0};
            int64_t lo = n * c / NCHUNK, hi = n * (c + 1) / NCHUNK;
            for (int64_t k = lo; k < hi; ++k) {
                if (!inside[k]) continue;
                const Sample& s = buf[k];
                int tap0 = (int)std::floor(s.mb) - 1;
                double ds = 0;
                for (int j = 0; j < 4; ++j) ds += bspline3_d(s.mb - (tap0 + j)) * Lmat[s.fb * B + tap0 + j];
                for (int j = 0; j < 3; ++j) {
                    g[j] += ds * s.dw[j];
                    g[3 + j] += ds * s.dt[j];
                }
            }
            std::copy(g, g + 6, part.data() + 6 * c);
        });
        for (int j = 0; j < 6; ++j) {
            double v = 0;
            for (int c = 0; c < NCHUNK; ++c) v += part[6 * c + j];
            grad_out[j] = v / total;
        }
        return total;
    }

    const Level& L;
    const Vol &F, &M;
    int B;
    double fscale, mscale;
    std::vector<double> hist, part, Lmat;
    std::vector<int64_t> count;
    std::vector<Sample> buf;
    std::vector<char> inside;
    Pool pool;
};

void skew_exp(const double* w, double* E) {
    double th = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    for (int i = 0; i < 9; ++i) E[i] = (i % 4 == 0);
    if (th < 1e-12) return;
    double k[3] = {w[0] / th, w[1] / th, w[2] / th};
    double K[9] = {0, -k[2], k[1], k[2], 0, -k[0], -k[1], k[0], 0};
    double KK[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) KK[3 * i + j] = K[3 * i] * K[j] + K[3 * i + 1] * K[3 + j] + K[3 * i + 2] * K[6 + j];
    double s = std::sin(th), c = 1 - std::cos(th);
    for (int i = 0; i < 9; ++i) E[i] += s * K[i] + c * KK[i];
}

// R <- exp([step[0:3]]x) R, t <- t + step[3:6]
void apply_step(const double* step, double* R, double* t) {
    double E[9], T[9];
    skew_exp(step, E);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) T[3 * i + j] = E[3 * i] * R[j] + E[3 * i + 1] * R[3 + j] + E[3 * i + 2] * R[6 + j];
    std::memcpy(R, T, sizeof T);
    for (int j = 0; j < 3; ++j) t[j] += step[3 + j];
}

}  // namespace

EXPORT int rigid_core_api_version() { return 1; }

// MI (and gradient, when grad is not null) on the given fixed indices u (n x 3).
EXPORT int64_t rigid_core_evaluate(const Level* L, const double* R, const double* t, const double* u, int64_t n,
                                   double* mi, double* grad) {
    Problem p(*L);
    return p.evaluate(R, t, 0, n, u, mi, grad);
}

// One pyramid level of adaptive stochastic gradient ascent. R (row-major 3x3) and t are updated in
// place; history receives (R, t) after every iteration, 12 doubles each.
// Returns 0, or -1 when too few samples map inside the moving image.
EXPORT int rigid_core_level(const Level* L, double* R, double* t, int iterations, int64_t samples,
                            double step_mm, double* history) {
    Problem p(*L);
    double rho = p.rho(0);
    double scale[6] = {1 / rho, 1 / rho, 1 / rho, 1, 1, 1};
    uint64_t draw = 1;

    auto gradient = [&](double* g) -> bool {
        double mi;
        for (int attempt = 0; attempt < 8; ++attempt)
            if (p.evaluate(R, t, draw++, samples, nullptr, &mi, g) >= 64) return true;
        return false;
    };

    // Gain: the first step moves points by about step_mm.
    double g[6], gnorm = 0;
    for (int k = 0; k < 5; ++k) {
        if (!gradient(g)) return -1;
        double s = 0;
        for (int j = 0; j < 6; ++j) s += (g[j] * scale[j]) * (g[j] * scale[j]);
        gnorm += std::sqrt(s) / 5;
    }
    double a = step_mm * std::pow(ASGD_A + 1, ASGD_ALPHA) / std::max(gnorm, 1e-12);
    double time = 0, prev[6];
    bool have_prev = false;
    for (int it = 0; it < iterations; ++it) {
        if (!gradient(g)) return -1;
        for (int j = 0; j < 6; ++j) g[j] *= scale[j];
        if (have_prev) {
            double dot = 0, n1 = 0, n2 = 0;
            for (int j = 0; j < 6; ++j) {
                dot += g[j] * prev[j];
                n1 += g[j] * g[j];
                n2 += prev[j] * prev[j];
            }
            double cosv = dot / (std::sqrt(n1) * std::sqrt(n2) + 1e-30);
            time = std::max(0.0, time + SIGMOID_MIN + (SIGMOID_MAX - SIGMOID_MIN) / (1 + std::exp(4 * cosv)));
        }
        std::memcpy(prev, g, sizeof g);
        have_prev = true;
        double gain = a / std::pow(ASGD_A + time + 1, ASGD_ALPHA);
        double gnorm2 = 0;
        for (int j = 0; j < 6; ++j) gnorm2 += g[j] * g[j];
        double moved = gain * std::sqrt(gnorm2);  // mm, rotation counted at the RMS radius
        if (moved > MAX_STEP_RATIO * step_mm) gain *= MAX_STEP_RATIO * step_mm / moved;
        double step[6];
        for (int j = 0; j < 6; ++j) step[j] = gain * g[j] * scale[j];
        apply_step(step, R, t);
        std::memcpy(history + 12 * it, R, 9 * sizeof(double));
        std::memcpy(history + 12 * it + 9, t, 3 * sizeof(double));
    }
    return 0;
}

// Deterministic finish: regular-step gradient ascent on the fixed sample set draw 0 (samples points).
// Returns the number of iterations taken, or -1 when too few samples map inside the moving image.
EXPORT int rigid_core_refine(const Level* L, double* R, double* t, int64_t samples, double step, double relax,
                             double min_step, int iterations) {
    Problem p(*L);
    double rho = p.rho(0);
    double scale[6] = {1 / rho, 1 / rho, 1 / rho, 1, 1, 1};
    double g[6], prev[6] = {0, 0, 0, 0, 0, 0}, mi;
    bool have_prev = false;
    int it = 0;
    for (; it < iterations; ++it) {
        if (p.evaluate(R, t, 0, samples, nullptr, &mi, g) < 64) return -1;
        double n = 0, dot = 0;
        for (int j = 0; j < 6; ++j) {
            g[j] *= scale[j];
            n += g[j] * g[j];
            dot += g[j] * prev[j];
        }
        if (have_prev && dot < 0) step *= relax;
        if (step < min_step) break;
        std::memcpy(prev, g, sizeof g);
        have_prev = true;
        double d[6];
        for (int j = 0; j < 6; ++j) d[j] = step * g[j] / (std::sqrt(n) + 1e-30) * scale[j];
        apply_step(d, R, t);
    }
    return it;
}
