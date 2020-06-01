//
// Created by lasagnaphil on 20. 6. 1..
//

#include "doctest.h"

#include "glmx/kinematic_map.h"
#include "glmx/common.h"

using namespace glm;
using namespace glmx;

TEST_CASE("euler <-> quat") {
    vec3 r = vec3(0.5, 1, 1.5);
    quat q = euler_to_quat(r);
    quat q2 = Rz(r.x) * Ry(r.y) * Rx(r.z);
    REQUIRE(length(log(inverse(q) * q2)) == doctest::Approx(0.0f));

    vec3 r2 = quat_to_euler(q);
    REQUIRE(r.x - r2.x == doctest::Approx(0.0f));
    REQUIRE(r.y - r2.y == doctest::Approx(0.0f));
    REQUIRE(r.z - r2.z == doctest::Approx(0.0f));
}


TEST_CASE("expmap jacobian") {
    vec3 r = vec3(1, 2, 3);
    vec3 dr = vec3(1e-3);
    float dt = 1e-3;
    quat q = exp(r);
    vec3 dq = log(inverse(q) * exp(r+dr));
    auto [J, Jdot] = expmap_body_jacobian_and_deriv(r, dr/dt);
    vec3 dq2 = J * dr;
    REQUIRE(dq.x - dq2.x == doctest::Approx(0.0f).epsilon(1e-4));
    REQUIRE(dq.y - dq2.y == doctest::Approx(0.0f).epsilon(1e-4));
    REQUIRE(dq.z - dq2.z == doctest::Approx(0.0f).epsilon(1e-4));

    auto Jp = expmap_body_jacobian(r+dr);
    mat3 Jdot2 = (Jp - J) / dt;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            REQUIRE(Jdot2[i][j] - Jdot[i][j] == doctest::Approx(0.0f).epsilon(1e-3));
        }
    }
}

TEST_CASE("euler jacobian") {
    vec3 r = vec3(0.5, 1, 1.5);
    vec3 dr = vec3(1e-3);
    float dt = 1.0 / 240.0f;
    quat q = euler_to_quat(r);
    vec3 dq = log(inverse(q) * euler_to_quat(r+dr));
    auto [J, Jdot] = euler_body_jacobian_and_deriv(r, dr/dt);
    vec3 dq2 = J * dr;
    REQUIRE(dq.x - dq2.x == doctest::Approx(0.0f).epsilon(1e-4));
    REQUIRE(dq.y - dq2.y == doctest::Approx(0.0f).epsilon(1e-4));
    REQUIRE(dq.z - dq2.z == doctest::Approx(0.0f).epsilon(1e-4));

    auto Jp = euler_body_jacobian(r+dr);
    mat3 Jdot2 = (Jp - J) / dt;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            REQUIRE(Jdot2[i][j] - Jdot[i][j] == doctest::Approx(0.0f).epsilon(1e-3));
        }
    }
}

