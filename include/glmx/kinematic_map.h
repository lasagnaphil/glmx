//
// Created by lasagnaphil on 20. 5. 26..
//

#ifndef GLMX_KINEMATIC_MAP_H
#define GLMX_KINEMATIC_MAP_H

// Calculates derivates of the exponential map.
// For more details, see these papers:
// - "Kinematic and dynamic modeling of spherical joints using exponential coordinates"
// - http://ethaneade.org/exp_diff.pdf

#include <glmx/common.h>
#include <glmx/spatial.h>

#include <tuple>

#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/trigonometric.hpp>

namespace glmx {
    inline std::tuple<float, float, float> expmap_deriv_params(float theta) {
        float alpha, beta, gamma;
        if (theta < 1e-4f) {
            alpha = 1 - theta*theta/6;
            beta = 0.5f - theta*theta/24;
            gamma = 1.f/6.f - theta*theta/120;
        }
        else {
            alpha = glm::sin(theta) / theta;
            beta = (1.0f - glm::cos(theta)) / (theta*theta);
            gamma = (1.0f - alpha) / (theta*theta);
        }
        return {alpha, beta, gamma};
    }

    inline std::tuple<float, float, float, float, float, float>
        expmap_deriv_params2(float theta, float d_theta) {

        float alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot;
        if (theta < 1e-4f) {
            alpha = 1 - theta*theta/6;
            beta = 0.5f - theta*theta/24;
            gamma = 1.f/6.f - theta*theta/120;
            alpha_dot = (-1.f/3.f + 1.f/30.f * theta*theta) * d_theta;
            beta_dot = (-1.f/12.f + 1.f/180.f * theta*theta) * d_theta;
            gamma_dot = (-1.f/60.f + 1.f/1260.f * theta*theta) * d_theta;
        }
        else {
            alpha = glm::sin(theta) / theta;
            beta = (1.0f - glm::cos(theta)) / (theta*theta);
            gamma = (1.0f - alpha) / (theta*theta);
            alpha_dot = (gamma - beta) * d_theta;
            beta_dot = (alpha - 2*beta) / (theta*theta) * d_theta;
            gamma_dot = (beta - 3*gamma) / (theta*theta) * d_theta;
        }
        return {alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot};
    }

    // Calculates the body jacobian of the exponential map of SO(3).
    inline glm::mat3 expmap_body_jacobian(glm::vec3 w) {
        float theta = glm::length(w);
        auto [alpha, beta, gamma] = expmap_deriv_params(theta);
        return alpha*glm::mat3(1.0f) - beta*skew_symmetric(w) + gamma*cartesian_product(w, w);
    }

    inline std::tuple<glm::mat3, glm::mat3> expmap_body_jacobian_and_deriv(glm::vec3 w, glm::vec3 w_dot) {
        float theta = glm::length(w);
        float w_w_dot = glm::dot(w, w_dot);
        auto [alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot] =
            expmap_deriv_params2(theta, w_w_dot);

        glm::mat3 J = alpha*glm::mat3(1.0f) - beta*skew_symmetric(w) + gamma*cartesian_product(w, w);
        glm::mat3 Jdot = alpha_dot*glm::mat3(1.0f) - beta_dot*skew_symmetric(w) - beta*skew_symmetric(w_dot)
                         + gamma_dot*cartesian_product(w, w)
                         + gamma*(cartesian_product(w, w_dot) + cartesian_product(w_dot, w));

        return {J, Jdot};
    }

    inline glm::mat3 spatial_jacobian(glm::vec3 w) {
        float theta = glm::length(w);
        auto [alpha, beta, gamma] = expmap_deriv_params(theta);
        return alpha*glm::mat3(1) + beta*skew_symmetric(w) + gamma*cartesian_product(w, w);
    }

    // Calculates the body derivative of the exponential map of SO(3).
    inline glm::vec3 expmap_body_deriv(glm::vec3 w, glm::vec3 w_dot) {
        float theta = glm::length(w);
        auto [alpha, beta, gamma] = expmap_deriv_params(theta);
        glm::vec3 w1 = glm::cross(w_dot, w);
        return alpha*w_dot - beta*glm::cross(w, w_dot) + gamma*glm::dot(w, w_dot)*w;
    }

    inline glm::vec3 expmap_spatial_deriv(glm::vec3 w, glm::vec3 w_dot) {
        float theta = glm::length(w);
        auto [alpha, beta, gamma] = expmap_deriv_params(theta);
        glm::vec3 w1 = glm::cross(w_dot, w);
        return alpha*w_dot + beta*glm::cross(w, w_dot) + gamma*glm::dot(w, w_dot)*w;
    }

    // Calculates the body derivative of the exponential map of SE(3).
    inline screw expmap_body_deriv(screw V, screw V_dot) {
        auto& w = V.w, v = V.v, w_dot = V_dot.w, v_dot = V_dot.v;
        float theta = glm::length(w);
        float w_v = glm::dot(w, v);
        auto [alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot] =
        expmap_deriv_params2(theta, w_v);
        screw Vdiff;
        Vdiff.w = alpha*w_dot - beta*glm::cross(w, w_dot) + gamma*glm::dot(w, w_dot)*w;
        Vdiff.v = -beta*glm::cross(v, v_dot) + gamma*(glm::dot(v, v_dot)*w + glm::dot(w, v_dot)*v)
                  + alpha_dot*v_dot - beta_dot*glm::cross(w, v_dot) + gamma_dot*glm::dot(w, v_dot)*w;
        return Vdiff;
    }

    // Calculates the second body derivative of the exponential map of SO(3).
    inline glm::vec3 expmap_second_body_deriv(glm::vec3 w, glm::vec3 w_dot, glm::vec3 w_2dot) {
        float theta = glm::length(w);
        float w_w_dot = glm::dot(w, w_dot);
        auto [alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot] =
        expmap_deriv_params2(theta, w_w_dot);
        return alpha*w_2dot
               + (gamma*w_w_dot + alpha_dot)*w_dot
               + glm::cross(beta_dot*w_dot+beta*w_2dot, w)
               + (gamma_dot*w_w_dot + gamma*(glm::length2(w_dot) + glm::dot(w, w_2dot)))*w;
    }

    inline screw expmap_second_body_deriv(screw V, screw V_dot, screw V_2dot) {
        fprintf(stderr, "Unimplemented!\n");
        exit(EXIT_FAILURE);
        return screw();
    }

    inline void reparameterize_expmap(glm::vec3& r, glm::vec3& r_dot, glm::vec3& r_2dot) {
        float theta = glm::length(r);
        const float pi = glm::pi<float>();
        if (theta <= pi) { return; }
        float eta = 1 - 2*pi/theta;
        float theta_3 = theta*theta*theta;
        float theta_5 = theta_3*theta*theta;
        float r_r_dot = glm::dot(r, r_dot);
        glm::vec3 rp = eta * r;
        glm::vec3 rp_dot = eta * r_dot + 2*pi*r_r_dot/theta_3*r;
        glm::vec3 rp_2dot = eta * r_2dot + 4*pi*r_r_dot/theta_3*r_dot
                            + 2*pi*((glm::length2(r_dot) + glm::dot(r, r_2dot))/theta_3 - 3*r_r_dot*r_r_dot/theta_5)*r;
        // glm::vec3 taup = tau - 2*pi/(eta*theta_3)*(glm::dot(r, tau)*r - glm::length2(r)*tau);
        r = rp; r_dot = rp_dot; r_2dot = rp_2dot;
    }

    // Calculates derivates of the ZYX euler map.

    inline glm::quat euler_to_quat(glm::vec3 q) {
        float s1 = glm::sin(0.5f*q.x);
        float s2 = glm::sin(0.5f*q.y);
        float s3 = glm::sin(0.5f*q.z);
        float c1 = glm::cos(0.5f*q.x);
        float c2 = glm::cos(0.5f*q.y);
        float c3 = glm::cos(0.5f*q.z);
        return glm::quat(c3*c2*c1 + s3*s2*s1, s3*c2*c1 - c3*s2*s1, c3*s2*c1 + s3*c2*s1, c3*c2*s1 - s3*s2*c1);
    }

    inline glm::vec3 quat_to_euler(glm::quat q) {
        glm::vec3 e;
        float sinr_cosp = 2*(q.w*q.x + q.y*q.z);
        float cosr_cosp = 1 - 2*(q.x*q.x+q.y*q.y);
        e.z = glm::atan(sinr_cosp, cosr_cosp);

        float sinp = 2*(q.w*q.y - q.z*q.x);
        if (glm::abs(sinp) >= 1) {
            e.y = M_PI/2 * glm::sign(sinp);
        }
        else {
            e.y = glm::asin(sinp);
        }

        float siny_cosp = 2*(q.w*q.z + q.x*q.y);
        float cosy_cosp = 1 - 2*(q.y*q.y + q.z*q.z);
        e.x = glm::atan(siny_cosp, cosy_cosp);
        return e;
    }

    inline glm::vec3 euler_to_expmap(glm::vec3 angles) {
        return log(euler_to_quat(angles));
    }

    inline glm::mat3 euler_body_jacobian(glm::vec3 q) {
        float s2 = glm::sin(q.y);
        float s3 = glm::sin(q.z);
        float c2 = glm::cos(q.y);
        float c3 = glm::cos(q.z);
        return glm::mat3(-s2, c2*s3, c2*c3, 0, c3, -s3, 1, 0, 0);
    }

    inline std::tuple<glm::mat3, glm::mat3> euler_body_jacobian_and_deriv(glm::vec3 q, glm::vec3 qdot) {
        float s2 = glm::sin(q.y);
        float s3 = glm::sin(q.z);
        float c2 = glm::cos(q.y);
        float c3 = glm::cos(q.z);
        glm::mat3 J = glm::mat3(-s2, c2*s3, c2*c3, 0, c3, -s3, 1, 0, 0);
        glm::mat3 Jdot = glm::mat3(-c2*qdot.y, -s2*s3*qdot.y + c2*c3*qdot.z, -s2*c3*qdot.y - c2*s3*qdot.z,
                                   0, -s3*qdot.z, -c3*qdot.z,
                                   0, 0, 0);
        return {J, Jdot};
    }

    inline void reparameterize_euler(glm::vec3& r) {
        const float pi = glm::pi<float>();
        if (r.x >= pi || r.x <= -pi) {
            int n = (int)glm::floor(glm::floor(r.x/pi + 1) / 2);
            r.x -= 2*n*pi;
        }
        if (r.y >= pi || r.y <= -pi) {
            int n = (int)glm::floor(glm::floor(r.y/pi + 1) / 2);
            r.y -= 2*n*pi;
        }
        if (r.z >= pi || r.z <= -pi) {
            int n = (int)glm::floor(glm::floor(r.z/pi + 1) / 2);
            r.z -= 2*n*pi;
        }
    }

    inline glm::vec3 rot_tracking_error(glm::vec3 r, glm::vec3 r_d) {
        float theta = length(r);
        float theta_d = length(r_d);
        float r_r_d = glm::dot(r, r_d);
        float alpha, beta, alpha_dot, beta_dot, alpha_d, beta_d;
        if (theta < 1e-4f) {
            alpha = 1 - theta*theta/6;
            beta = 0.5f - theta*theta/24;
            alpha_dot = -1.f/3.f + 1.f/30.f * theta*theta;
            beta_dot = -1.f/12.f + 1.f/180.f * theta*theta;
            alpha_d = 1 - theta_d*theta_d/6;
            beta_d = 0.5f - theta_d*theta_d/24;
        }
        else {
            alpha = glm::sin(theta) / theta;
            beta = (1.0f - glm::cos(theta)) / (theta*theta);
            float gamma = (1.0f - alpha) / (theta*theta);
            alpha_dot = gamma - beta;
            beta_dot = (alpha - 2*beta) / (theta*theta);
            alpha_d = glm::sin(theta_d) / theta_d;
            beta_d = (1.0f - glm::cos(theta_d)) / (theta_d*theta_d);
        }
        return (0.5f*(1.0f+glm::cos(theta_d))*alpha - alpha_dot*alpha_d*r_r_d - 0.5f*beta_d*beta_dot*r_r_d*r_r_d)*r
               - (alpha_d*alpha + beta_d*beta*r_r_d)*r_d;
    }
}

#endif //GLMX_KINEMATIC_MAP_H
