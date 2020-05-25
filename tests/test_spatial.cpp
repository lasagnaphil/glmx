//
// Created by lasagnaphil on 20. 5. 22..
//

#include "doctest.h"
#include "glmx/spatial.h"

using namespace glm;
using namespace glmx;

TEST_CASE("sym_mat3 * vec3") {
    sym_mat3 I(1, 2, 3, 4, 5, 6);
    vec3 v(1, 2, 3);
    vec3 v1 = I * v;
    vec3 v2 = mat3_cast(I) * v;
    REQUIRE(v1 == v2);
}

TEST_CASE("move_frame with sym_mat3") {
    sym_mat3 I_b(1, 2, 3, 4, 5, 6);
    mat3 I_b_ns = mat3_cast(I_b);
    quat q_ba = Rx(0.1) * Ry(0.2) * Rz(0.3);
    mat3 R_ba = mat3_cast(q_ba);
    sym_mat3 I_a = move_frame(I_b, q_ba);
    mat3 I_a_ns = transpose(R_ba) * I_b_ns * R_ba;
    REQUIRE(I_a.xx == doctest::Approx(I_a_ns[0][0]));
    REQUIRE(I_a.yy == doctest::Approx(I_a_ns[1][1]));
    REQUIRE(I_a.zz == doctest::Approx(I_a_ns[2][2]));
    REQUIRE(I_a.xy == doctest::Approx(I_a_ns[0][1]));
    REQUIRE(I_a.yz == doctest::Approx(I_a_ns[1][2]));
    REQUIRE(I_a.zx == doctest::Approx(I_a_ns[2][0]));
}

TEST_CASE("quadratic_form with sym_mat3") {
    sym_mat3 I(1, 2, 3, 4, 5, 6);
    vec3 w(1, 2, 3);
    vec3 v(4, 5, 6);
    REQUIRE(quadratic_form(I, w, v) == quadratic_form(mat3_cast(I), w, v));
    REQUIRE(quadratic_form(I, w) == quadratic_form(I, w, w));
    REQUIRE(quadratic_form(I, w) == quadratic_form(mat3_cast(I), w));
}

TEST_CASE("symmetric_cartesian_product with sym_mat3") {
    vec3 w(1, 2, 3);
    sym_mat3 I = symmetric_cartesian_product(w);
    mat3 I_ns = cartesian_product(w, w);
    REQUIRE(I.xx == I_ns[0][0]);
    REQUIRE(I.yy == I_ns[1][1]);
    REQUIRE(I.zz == I_ns[2][2]);
    REQUIRE(I.xy == I_ns[0][1]);
    REQUIRE(I.yz == I_ns[1][2]);
    REQUIRE(I.zx == I_ns[2][0]);
}

TEST_CASE("inverse with sym_mat3") {
    sym_mat3 I(1, 2, 3, 4, 5, 6);
    sym_mat3 I_inv = inverse(I);
    mat3 I_ns = mat3_cast(I);
    mat3 I_inv_ns = glm::inverse(I_ns);
    REQUIRE(I_inv.xx == I_inv_ns[0][0]);
    REQUIRE(I_inv.yy == I_inv_ns[1][1]);
    REQUIRE(I_inv.zz == I_inv_ns[2][2]);
    REQUIRE(I_inv.xy == I_inv_ns[0][1]);
    REQUIRE(I_inv.yz == I_inv_ns[1][2]);
    REQUIRE(I_inv.zx == I_inv_ns[2][0]);
}

TEST_CASE("move_frame with spatial_inertia_mat") {
    spatial_inertia_mat G_b(sym_mat3(1, 2, 3, 4, 5, 6), glm::vec3(7, 8, 9), 10);
    transform T_ba(vec3(1, 2, 3), Rx(0.1)*Ry(0.2)*Rz(0.3));
    spatial_inertia_mat G_a = move_frame(G_b, T_ba);
    sym_mat6 G_a_ns = move_frame(sym_mat6(G_b), T_ba);
    vec3 c_ns = skew_symmetric_cast(G_a_ns.C) / G_a_ns.M.xx;
    REQUIRE(G_a.I.xx == doctest::Approx(G_a_ns.I.xx));
    REQUIRE(G_a.I.yy == doctest::Approx(G_a_ns.I.yy));
    REQUIRE(G_a.I.zz == doctest::Approx(G_a_ns.I.zz));
    REQUIRE(G_a.I.xy == doctest::Approx(G_a_ns.I.xy));
    REQUIRE(G_a.I.yz == doctest::Approx(G_a_ns.I.yz));
    REQUIRE(G_a.I.zx == doctest::Approx(G_a_ns.I.zx));
    REQUIRE(G_a.c.x == doctest::Approx(c_ns.x));
    REQUIRE(G_a.c.y == doctest::Approx(c_ns.y));
    REQUIRE(G_a.c.z == doctest::Approx(c_ns.z));
    REQUIRE(G_a.m == doctest::Approx(G_a_ns.M.xx));
}
