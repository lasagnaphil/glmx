//
// Created by lasagnaphil on 19. 11. 7..
//

#ifndef GENGINE_QUAT_H
#define GENGINE_QUAT_H

#include <glm/vec3.hpp>
#include <glm/gtx/quaternion.hpp>
#include "euler.h"

namespace glmx {
    inline glm::vec3 log(glm::quat q) {
        constexpr float pi = glm::pi<float>();
        q = glm::normalize(q);
        float a = glm::sqrt(1 - q.w*q.w);
        if (a <= glm::epsilon<float>()) {
            return glm::vec3(0);
        }
        float theta = 2.0f * glm::atan(a, q.w);
        if (theta > pi) {
            theta -= 2*pi;
        }
        else if (theta < -pi) {
            theta += 2*pi;
        }
        glm::vec3 v = theta / a * glm::vec3(q.x, q.y, q.z);
        return v;
    }

    inline glm::vec3 logdiff(glm::quat q1, glm::quat q2) {
        // return glmx::log(q2 * glm::conjugate(q1));
        return glmx::log(glm::conjugate(q1) * q2);
    }

    inline glm::quat exp(glm::vec3 v) {
        float theta = glm::length(v);
        if (theta <= glm::epsilon<float>()) {
            return glm::identity<glm::quat>();
        }
        glm::vec3 u = v / theta;
        return glm::quat(glm::cos(theta/2), glm::sin(theta/2) * u);
    }

    inline float extractXRot(glm::quat q) {
        if (q.x * q.x + q.w * q.w <= glm::epsilon<float>()) {
            return 0.0f;
        }
        return 2 * glm::atan(q.x, q.w);
    }

    inline float extractYRot(glm::quat q) {
        if (q.y * q.y + q.w * q.w <= glm::epsilon<float>()) {
            return 0.0f;
        }
        return 2 * glm::atan(q.y, q.w);
    }

    inline float extractZRot(glm::quat q) {
        if (q.z * q.z + q.w * q.w <= glm::epsilon<float>()) {
            return 0.0f;
        }
        return 2 * glm::atan(q.z, q.w);
    }

    inline glm::mat4 rotMatrixBetweenVecs(glm::vec3 a, glm::vec3 b) {
        glm::vec3 v = glm::cross(a, b);
        float s2 = glm::dot(v, v);
        if (s2 < glm::epsilon<float>()) {
            return glm::mat4(1.0f);
        }
        else {
            // Rodrigue's formula
            float c = glm::dot(a, b);
            glm::mat3 vhat;
            vhat[0][0] = vhat[1][1] = vhat[2][2] = 0;
            vhat[2][1] = v[0]; vhat[1][2] = -v[0];
            vhat[0][2] = v[1]; vhat[2][0] = -v[1];
            vhat[1][0] = v[2]; vhat[0][1] = -v[2];
            return glm::mat3(1.0f) + vhat + vhat*vhat*(1 - c)/(s2);
        }
    }

    // https://math.stackexchange.com/questions/90081/quaternion-distance
    inline glm::quat quatBetweenVecs(glm::vec3 a, glm::vec3 b) {
        glm::vec3 w = glm::cross(a, b);
        glm::quat q = glm::quat(1.f + dot(a, b), w);
        return glm::normalize(q);
    }

    inline float angleBetweenQuats(glm::quat q1, glm::quat q2) {
        float inner = glm::dot(q1, q2);
        float angle = glm::acos(2*inner*inner - 1);
        return inner;
    }

    inline glm::vec3 Ex() { return glm::vec3(1, 0, 0); }
    inline glm::vec3 Ey() { return glm::vec3(0, 1, 0); }
    inline glm::vec3 Ez() { return glm::vec3(0, 0, 1); }

    inline glm::quat Rx(float theta) { return glm::angleAxis(theta, Ex()); }
    inline glm::quat Ry(float theta) { return glm::angleAxis(theta, Ey()); }
    inline glm::quat Rz(float theta) { return glm::angleAxis(theta, Ez()); }

    inline glm::mat3 mat3_from_diag(glm::vec3 v) {
        return glm::mat3(v.x, 0, 0, 0, v.y, 0, 0, 0, v.z);
    }

    inline glm::mat3 cartesian_product(glm::vec3 w, glm::vec3 v) {
        glm::mat3 m;
        m[0][0] = w[0]*v[0];
        m[0][1] = w[1]*v[0];
        m[0][2] = w[2]*v[0];
        m[1][0] = w[0]*v[1];
        m[1][1] = w[1]*v[1];
        m[1][2] = w[2]*v[1];
        m[2][0] = w[0]*v[2];
        m[2][1] = w[1]*v[2];
        m[2][2] = w[2]*v[2];
        return m;
    }

    inline glm::mat3 skew_symmetric(glm::vec3 w) {
        glm::mat3 m(0.0f);
        m[1][2] = w.x;
        m[2][1] = -w.x;
        m[2][0] = w.y;
        m[0][2] = -w.y;
        m[0][1] = w.z;
        m[1][0] = -w.z;
        return m;
    }

    inline glm::vec3 skew_symmetric_cast(glm::mat3 m) {
        return glm::vec3(m[1][2], m[2][0], m[0][1]);
    }
}

#endif //GENGINE_QUAT_H
