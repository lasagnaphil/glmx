//
// Created by lasagnaphil on 20. 4. 29..
//

#ifndef ARTSIM_MATH_H
#define ARTSIM_MATH_H

#include <glm/vec3.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <tuple>

#include "glmx/transform.h"
#include "glmx/common.h"

namespace glmx {
    struct screw {
        glm::vec3 w;
        glm::vec3 v;

        screw() : w(0.0f), v(0.0f) {}
        screw(glm::vec3 w, glm::vec3 v) : w(w), v(v) {}
    };

    inline screw make_screw(const float* ptr) {
        return screw(glm::vec3(ptr[0], ptr[1], ptr[2]), glm::vec3(ptr[3], ptr[4], ptr[5]));
    }

    inline screw operator+(const screw& V1, const screw& V2) {
        return {V1.w + V2.w, V1.v + V2.v};
    }

    inline screw& operator+=(screw& V1, const screw& V2) {
        V1.w += V2.w; V1.v += V2.v;
        return V1;
    }

    inline screw operator-(const screw& V1, const screw& V2) {
        return {V1.w - V2.w, V1.v - V2.v};
    }

    inline screw& operator-=(screw& V1, const screw& V2) {
        V1.w -= V2.w; V1.v -= V2.v;
        return V1;
    }

    inline screw operator*(const screw& V, float theta) {
        return {theta * V.w, theta * V.v};
    }

    inline screw operator*(float theta, const screw& V) {
        return {theta * V.w, theta * V.v};
    }

    inline screw operator/(const screw& V, float theta) {
        return {V.w / theta, V.v / theta};
    }

    inline screw operator-(const screw& V) {
        return {-V.w, -V.v};
    }

    inline float dot(const screw& V1, const screw& V2) {
        return glm::dot(V1.w, V2.w) + glm::dot(V1.v, V2.v);
    }

    inline float length(const screw& V) {
        float w2 = glm::length2(V.w);
        if (w2 >= glm::epsilon<float>()) {
            return glm::sqrt(w2);
        }
        float v2 = glm::length2(V.v);
        if (v2 >= glm::epsilon<float>()) {
            return glm::sqrt(v2);
        }
        return 0.0f;
    }

    inline screw normalize(const screw& V) {
        float w2 = glm::length2(V.w);
        if (w2 >= glm::epsilon<float>()) {
            return V / glm::sqrt(w2);
        }
        float v2 = glm::length2(V.v);
        if (v2 >= glm::epsilon<float>()) {
            return V / glm::sqrt(v2);
        }
        return screw();
    }

    // Calculates exp([V] * theta).
    inline transform screw_move(screw V, float theta) {
        glm::vec3 w_cross_v = glm::cross(V.w, V.v);
        glm::vec3 p = V.v * theta + (1 - glm::cos(theta)) * w_cross_v +
                      (theta - glm::sin(theta)) * glm::cross(V.w, w_cross_v);

        return transform(p, exp(V.w * theta));
    }

    inline transform exp(const screw& V) {
        float V_len = length(V);
        if (V_len < glm::epsilon<float>()) {
            return transform();
        }
        screw V_hat = V / V_len;
        return screw_move(V_hat, V_len);
    }

    inline screw log(const transform& T) {
        screw V;
        V.w = log(T.q);
        float theta = glm::length(V.w);
        if (theta <= glm::epsilon<float>()) {
            V.w = glm::vec3(0);
            V.v = T.v;
        }
        else {
            float delta;
            if (theta < 1e-4f) {
                delta = 12.0f;
            }
            else {
                float alpha = sin(theta) / theta;
                float beta = (1 - cos(theta)) / (theta*theta);
                delta = (1 - alpha/(2.0f*beta)) / (theta*theta);
            }
            glm::vec3 w_cross_v = glm::cross(V.w, T.v);
            V.v = T.v - 0.5f * w_cross_v + delta*glm::cross(V.w, w_cross_v);
        }
        return V;
    }

    /*
    inline transform adjoint_to_trans(transform T, screw V) {
        glm::vec3 w = T.q * V.w;
        return transform(glm::cross(T.v, w) + T.q * V.v, artsim::exp(w));
    }
     */

    inline screw big_adj(transform T, screw V) {
        glm::vec3 w = T.q * V.w;
        return screw(w, glm::cross(T.v, w) + T.q * V.v);
    }

    inline screw big_T_adj(transform T, screw V) {
        return screw(glm::conjugate(T.q) * (V.w + glm::cross(V.v, T.v)), glm::conjugate(T.q) * V.v);
    }

    inline screw small_adj(screw V1, screw V2) {
        return screw(glm::cross(V1.w, V2.w), glm::cross(V1.v, V2.w) + glm::cross(V1.w, V2.v));
    }

    inline screw small_T_adj(screw V1, screw V2) {
        return screw(glm::cross(V2.w, V1.w) + glm::cross(V2.v, V1.v), glm::cross(V2.v, V1.w));
    }

    struct sym_mat3 {
        float xx;
        float yy;
        float zz;
        float xy;
        float yz;
        float zx;

        sym_mat3() = default;
        sym_mat3(float xx, float yy, float zz, float xy, float yz, float zx)
                : xx(xx), yy(yy), zz(zz), xy(xy), yz(yz), zx(zx) {}
        explicit sym_mat3(const glm::mat3& m)
                : xx(m[0][0]), yy(m[1][1]), zz(m[2][2]), xy(m[0][1]), yz(m[1][2]), zx(m[2][0]) {}
    };

    inline sym_mat3 operator+(const sym_mat3& I1, const sym_mat3& I2) {
        return sym_mat3(I1.xx + I2.xx, I1.yy + I2.yy, I1.zz + I2.zz, I1.xy + I2.xy, I1.yz + I2.yz, I1.zx + I2.zx);
    }

    inline sym_mat3& operator+=(sym_mat3& I1, const sym_mat3& I2) {
        I1.xx += I2.xx; I1.yy += I2.yy; I1.zz += I2.zz;
        I1.xy += I2.xy; I1.yz += I2.yz; I1.zx += I2.zx;
        return I1;
    }

    inline sym_mat3 operator-(const sym_mat3& I1, const sym_mat3& I2) {
        return sym_mat3(I1.xx - I2.xx, I1.yy - I2.yy, I1.zz - I2.zz, I1.xy - I2.xy, I1.yz - I2.yz, I1.zx - I2.zx);
    }

    inline sym_mat3& operator-=(sym_mat3& I1, const sym_mat3& I2) {
        I1.xx -= I2.xx; I1.yy -= I2.yy; I1.zz -= I2.zz;
        I1.xy -= I2.xy; I1.yz -= I2.yz; I1.zx -= I2.zx;
        return I1;
    }

    inline glm::vec3 operator*(const sym_mat3& I, glm::vec3 w) {
        return glm::vec3(
                I.xx*w.x + I.xy*w.y + I.zx*w.z,
                I.xy*w.x + I.yy*w.y + I.yz*w.z,
                I.zx*w.x + I.yz*w.y + I.zz*w.z);
    }

    inline sym_mat3 operator*(float a, const sym_mat3& I) {
        return sym_mat3(a*I.xx, a*I.yy, a*I.zz, a*I.xy, a*I.yz, a*I.zx);
    }

    inline glm::mat3 mat3_cast(const sym_mat3& I) {
        return glm::mat3(I.xx, I.xy, I.zx, I.xy, I.yy, I.yz, I.zx, I.yz, I.zz);
    }

    inline glm::mat3 operator*(const sym_mat3& I1, const sym_mat3& I2) {
        // TODO: This could be implemented more efficiently
        return mat3_cast(I1) * mat3_cast(I2);
    }

    inline sym_mat3 move_frame(sym_mat3 I_b, glm::quat q_ba) {
        // TODO: This could be implemented more efficiently
        glm::mat3 R = glm::mat3_cast(q_ba);
        return sym_mat3(glm::transpose(R) * mat3_cast(I_b) * R);
    }

    inline float quadratic_form(const sym_mat3& I, glm::vec3 w) {
        return I.xx*w.x*w.x + I.yy*w.y*w.y + I.zz*w.z*w.z
            + 2*I.xy*w.x*w.y + 2*I.yz*w.y*w.z + 2*I.zx*w.z*w.x;
    }

    inline float quadratic_form(const sym_mat3& I, glm::vec3 w, glm::vec3 v) {
        return I.xx*w.x*v.x + I.yy*w.y*v.y + I.zz*w.z*v.z
            + I.xy*(w.x*v.y+w.y*v.x) + I.yz*(w.y*v.z+w.z*v.y) + I.zx*(w.z*v.x+w.x*v.z);
    }

    inline float quadratic_form(const glm::mat3& I, glm::vec3 w, glm::vec3 v) {
        return I[0][0]*w[0]*v[0] + I[0][1]*w[1]*v[0] + I[0][2]*w[2]*v[0]
             + I[1][0]*w[0]*v[1] + I[1][1]*w[1]*v[1] + I[1][2]*w[2]*v[1]
             + I[2][0]*w[0]*v[2] + I[2][1]*w[1]*v[2] + I[2][2]*w[2]*v[2];
    }

    inline float quadratic_form(const glm::mat3& I, glm::vec3 w) {
        return I[0][0]*w[0]*w[0] + I[0][1]*w[1]*w[0] + I[0][2]*w[2]*w[0]
             + I[1][0]*w[0]*w[1] + I[1][1]*w[1]*w[1] + I[1][2]*w[2]*w[1]
             + I[2][0]*w[0]*w[2] + I[2][1]*w[1]*w[2] + I[2][2]*w[2]*w[2];
    }

    inline sym_mat3 symmetric_cartesian_product(glm::vec3 w) {
        return sym_mat3(w.x * w.x, w.y * w.y, w.z * w.z, w.x * w.y, w.y * w.z, w.z * w.x);
    }

    inline float determinant(const sym_mat3& I) {
        return I.xx*I.yy*I.zz + 2*I.xy*I.yz*I.zx - I.xx*I.yz*I.yz - I.yy*I.zx*I.zx - I.zz*I.xy*I.xy;
    }
    inline sym_mat3 cofactor(const sym_mat3& I) {
        return sym_mat3(I.yy*I.zz-I.yz*I.yz, I.zz*I.xx-I.zx*I.zx, I.xx*I.yy-I.xy*I.xy,
                        I.yz*I.zx-I.xy*I.zz, I.xy*I.zx-I.xx*I.yz, I.xy*I.yz-I.yy*I.zx);
    }

    inline sym_mat3 inverse(const sym_mat3& I) {
        return (1.0f / determinant(I)) * cofactor(I);
    }
    // Spatial matrix.
    /*
     *      Spatial mass matrix
           ---------- ----------
        0 |          |          |
        1 |  I(=mG)  |   m c×   |
        2 |          |          |
           ---------- ----------
        3 |          |          |
        4 |  -m c×   |    mE    |
        5 |          |          |
           ---------- ----------
     */
    struct spatial_inertia_mat {
        sym_mat3 I;
        glm::vec3 c;
        float m;

        spatial_inertia_mat() = default;
        spatial_inertia_mat(sym_mat3 I, glm::vec3 c, float m) : I(I), c(c), m(m) {}
    };

    inline spatial_inertia_mat operator+(const spatial_inertia_mat& G1, const spatial_inertia_mat& G2) {
        return {G1.I + G2.I, G1.c + G2.c, G1.m + G2.m};
    }

    inline spatial_inertia_mat& operator+=(spatial_inertia_mat& G1, const spatial_inertia_mat& G2) {
        G1.I += G2.I; G1.c += G2.c; G1.m += G2.m;
        return G1;
    }

    inline spatial_inertia_mat operator-(const spatial_inertia_mat& G1, const spatial_inertia_mat& G2) {
        return {G1.I - G2.I, G1.c - G2.c, G1.m - G2.m};
    }

    inline spatial_inertia_mat& operator-=(spatial_inertia_mat& G1, const spatial_inertia_mat& G2) {
        G1.I -= G2.I; G1.c -= G2.c; G1.m -= G2.m;
        return G1;
    }

    inline screw operator*(const spatial_inertia_mat& G, screw V) {
        return screw(G.I * V.w + G.m * glm::cross(G.c, V.v), G.m*(V.v - glm::cross(G.c, V.w)));
    }

    inline spatial_inertia_mat move_frame(const spatial_inertia_mat& G_b, transform T_ba) {
        spatial_inertia_mat G_a;
        glm::vec3 c = glm::conjugate(T_ba.q) * G_b.c;
        glm::vec3 cp = glm::conjugate(T_ba.q) * (G_b.c - T_ba.v);
        G_a.I = move_frame(G_b.I, T_ba.q);
        G_a.I.xx += G_b.m*(-c.y*c.y - c.z*c.z + cp.y*cp.y + cp.z*cp.z);
        G_a.I.yy += G_b.m*(-c.z*c.z - c.x*c.x + cp.z*cp.z + cp.x*cp.x);
        G_a.I.zz += G_b.m*(-c.x*c.x - c.y*c.y + cp.x*cp.x + cp.y*cp.y);
        G_a.I.xy += G_b.m*(c.x*c.y - cp.x*cp.y);
        G_a.I.yz += G_b.m*(c.y*c.z - cp.y*cp.z);
        G_a.I.zx += G_b.m*(c.z*c.x - cp.z*cp.x);
        G_a.c = cp;
        G_a.m = G_b.m;
        return G_a;
    }

    inline float quadratic_form(const spatial_inertia_mat& G, screw V) {
        return quadratic_form(G.I, V.w) + G.m*glm::length2(V.v) - 2*G.m*glm::dot(glm::cross(V.w, V.v), G.c);
    }

    struct sym_mat6 {
        sym_mat3 I, M;
        glm::mat3 C;

        sym_mat6() = default;
        sym_mat6(sym_mat3 I, sym_mat3 M, glm::mat3 C) : I(I), M(M), C(C) {}
        explicit sym_mat6(const spatial_inertia_mat& G)
            : I(G.I), M(sym_mat3(G.m, G.m, G.m, 0, 0, 0)), C(G.m * skew_symmetric(G.c)) {}
    };

    inline sym_mat6 operator+(const sym_mat6& G1, const sym_mat6& G2) {
        return sym_mat6(G1.I + G2.I, G1.M + G2.M, G1.C + G2.C);
    }

    inline sym_mat6& operator+=(sym_mat6& G1, const sym_mat6& G2) {
        G1.I += G2.I; G1.M += G2.M; G1.C += G2.C;
        return G1;
    }

    inline sym_mat6 operator-(const sym_mat6& G1, const sym_mat6& G2) {
        return sym_mat6(G1.I - G2.I, G1.M - G2.M, G1.C - G2.C);
    }

    inline sym_mat6& operator-=(sym_mat6& G1, const sym_mat6& G2) {
        G1.I -= G2.I; G1.M -= G2.M; G1.C -= G2.C;
        return G1;
    }

    inline screw operator*(const sym_mat6& G, screw V) {
        return screw(G.I * V.w + G.C * V.v, glm::transpose(G.C) * V.w + G.M * V.v);
    }

    inline sym_mat6 operator*(float a, const sym_mat6& G) {
        return sym_mat6(a*G.I, a*G.M, a*G.C);
    }

    inline sym_mat6 symmetric_cartesian_product(screw V) {
        return sym_mat6(symmetric_cartesian_product(V.w),
                        symmetric_cartesian_product(V.v),
                        cartesian_product(V.w, V.v));
    }


    inline sym_mat6 move_frame(const sym_mat6& G_b, transform T_ba) {
        sym_mat6 G_a;
        glm::mat3 P = skew_symmetric(T_ba.v);
        glm::mat3 R = glm::mat3_cast(T_ba.q);
        glm::mat3 I = mat3_cast(G_b.I);
        glm::mat3 PM = P*mat3_cast(G_b.M);
        glm::mat3 CP = G_b.C * P;
        G_a.I = move_frame(sym_mat3(I + CP + glm::transpose(CP) - PM*P), R);
        G_a.C = glm::transpose(R)*(G_b.C - PM)*R;
        G_a.M = move_frame(G_b.M, R);
        return G_a;
    }

    inline float quadratic_form(const sym_mat6& G, screw V) {
        return quadratic_form(G.I, V.w) + quadratic_form(G.M, V.v) + 2*quadratic_form(G.C, V.w, V.v);
    }

    inline sym_mat6 inverse(const sym_mat6& G) {
        glm::mat3 Iinv = mat3_cast(inverse(G.I));
        glm::mat3 Iinv_C = Iinv * G.C;
        glm::mat3 D = inverse(mat3_cast(G.M) - glm::transpose(G.C) * Iinv_C);
        sym_mat6 Ginv;
        Ginv.I = sym_mat3(Iinv + Iinv_C * D * glm::transpose(Iinv_C));
        Ginv.C = -Iinv_C * D;
        Ginv.M = sym_mat3(D);
        return Ginv;
    }
}


#endif //ARTSIM_MATH_H
