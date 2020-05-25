//
// Created by lasagnaphil on 20. 5. 26..
//

#ifndef GLMX_KINEMATIC_MAP_H
#define GLMX_KINEMATIC_MAP_H

// Calculates derivates of the exponential map.
// For more details, see these papers:
// - "Kinematic and dynamic modeling of spherical joints using exponential coordinates"
// - http://ethaneade.org/exp_diff.pdf

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
        return alpha*glm::mat3(1) - beta*skew_symmetric(w) + gamma*cartesian_product(w, w);
    }

    inline std::tuple<glm::mat3, glm::mat3> expmap_body_jacobian_and_deriv(glm::vec3 w, glm::vec3 w_dot) {
        float theta = glm::length(w);
        float w_w_dot = glm::dot(w, w_dot);
        auto [alpha, beta, gamma, alpha_dot, beta_dot, gamma_dot] =
        expmap_deriv_params2(theta, w_w_dot);

        glm::mat3 J = alpha*glm::mat3(1) - beta*skew_symmetric(w) + gamma*cartesian_product(w, w);
        glm::mat3 Jdot = -beta_dot*skew_symmetric(w) - beta*skew_symmetric(w_dot)
                         + gamma_dot*cartesian_product(w, w)
                         + (gamma - beta)*w_w_dot*glm::mat3(1)
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
        int n = -(int)glm::floor(glm::floor(theta/pi + 1) / 2);
        float eta = (1 + 2*n*pi/theta);
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
    // For more details, see the RBDL paper.

    inline glm::quat euler_to_quat(glm::vec3 q) {
        float s1 = glm::sin(0.5f*q.x);
        float s2 = glm::sin(0.5f*q.y);
        float s3 = glm::sin(0.5f*q.z);
        float c1 = glm::cos(0.5f*q.x);
        float c2 = glm::cos(0.5f*q.y);
        float c3 = glm::cos(0.5f*q.z);
        return glm::quat(c1*c2*c3 - s1*s2*s3, c1*c2*s3+s1*s2*c3, c1*s2*c3-s1*c2*s3, c1*s2*s3+s1*c2*c3);
    }

    inline glm::vec3 euler_to_expmap(glm::vec3 angles) {
        return log(euler_to_quat(angles));
    }

    inline std::tuple<glm::mat3, glm::mat3> euler_body_jacobian_and_deriv(glm::vec3 q, glm::vec3 qdot) {
        float s2 = glm::sin(q.y);
        float s3 = glm::sin(q.z);
        float c2 = glm::cos(q.y);
        float c3 = glm::cos(q.z);
        glm::mat3 J = glm::mat3(-s2, c2*s3, c3*c3, 0, c3, -s3, 1, 0, 0);
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
}

#endif //GLMX_KINEMATIC_MAP_H
