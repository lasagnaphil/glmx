//
// Created by lasagnaphil on 19. 10. 1..
//

#ifndef GENGINE_GLMX_TRANSFORM_H
#define GENGINE_GLMX_TRANSFORM_H

#include <glm/vec3.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>

namespace glmx {
    struct transform {
        glm::vec3 v;
        glm::quat q;

        transform() : v(0.0f), q(glm::identity<glm::quat>()) {}
        transform(glm::vec3 v) : v(v), q(glm::identity<glm::quat>()) {}
        transform(glm::quat q) : v(glm::vec3()), q(q) {}
        transform(glm::vec3 v, glm::quat q) : v(v), q(q) {}
    };

    inline glm::mat4 mat4_cast(const transform &t) {
        return glm::translate(t.v) * mat4_cast(t.q);
    }

    inline transform operator*(const transform &t1, const transform &t2) {
        return {t1.q * t2.v + t1.v, t1.q * t2.q};
    }

    inline transform operator*(glm::quat q, const transform& t) {
        return {q * t.v, q * t.q};
    }

    inline transform operator*(glm::vec3 v, const transform& t) {
        return {t.v + v, t.q};
    }

    inline transform operator/(const transform &t1, const transform &t2) {
        return {conjugate(t2.q) * (t1.v - t2.v), conjugate(t2.q) * t1.q};
    }

    inline transform inverse(const transform& t) {
        return transform(conjugate(t.q) * (-t.v), conjugate(t.q));
    }
}

#endif //GENGINE_TRANSFORM_H
