#include <iostream>
#include "RedMath.h"

int main()
{
    std::cout << "Hello there!\n";

    //-------------------------------------------------------------------   Vec4
    Vec4<float> vec4f_0{};
    Vec4<float> vec4f_1{15.0f};
    Vec4<float> vec4f_2{-10.0f, 3.14f};
    Vec4<float> vec4f_3{1.0f, 2.0f, 3.0f};
    Vec4<float> vec4f_4{33.0f, -500.0f, 200.0f, -0.1f};

    vec4f_0 = vec4f_1 + vec4f_2;
    vec4f_0 += vec4f_1 + vec4f_2;

    vec4f_0 = vec4f_1 - vec4f_2;
    vec4f_0 -= vec4f_1 - vec4f_2;

    vec4f_0 = vec4f_1 * vec4f_2;
    vec4f_0 *= vec4f_1 * vec4f_2;

    vec4f_0 = vec4f_1 / vec4f_4;
    vec4f_0 /= vec4f_1 / vec4f_4;

    float dot = vec4f_0.dot(vec4f_1);
    float length2 = vec4f_0.length2();
    float length = vec4f_0.length();

    Vec4<float> vec4f_5 = vec4f_0.normalized();
    vec4f_0.normalize();
    //-------------------------------------------------------------------   Vec4



    //-------------------------------------------------------------------   Vec3
    Vec3<float> vec3f_0{};
    Vec3<float> vec3f_1{ 15.0f };
    Vec3<float> vec3f_2{ -10.0f, 3.14f};
    Vec3<float> vec3f_3{ 1.0f, 2.0f, 3.0f };

    vec3f_0 = vec3f_1 + vec3f_2;
    vec3f_0 += vec3f_1 + vec3f_2;

    vec3f_0 = vec3f_1 - vec3f_2;
    vec3f_0 -= vec3f_1 - vec3f_2;

    vec3f_0 = vec3f_1 * vec3f_2;
    vec3f_0 *= vec3f_1 * vec3f_2;

    vec3f_0 = vec3f_1 / vec3f_2;
    vec3f_0 /= vec3f_1 / vec3f_2;

    dot = vec3f_0.dot(vec3f_1);
    length2 = vec3f_0.length2();
    length = vec3f_0.length();

    Vec3<float> vec3f_4 = vec3f_0.normalized();
    vec3f_0.normalize();
    //-------------------------------------------------------------------   Vec3



	//-------------------------------------------------------------------   Vec2
    Vec2<float> vec2f_0{};
    Vec2<float> vec2f_1{ 15.0f };
    Vec2<float> vec2f_2{ -10.0f, 3.14f };

    vec2f_0 = vec2f_1 + vec2f_2;
    vec2f_0 += vec2f_1 + vec2f_2;

    vec2f_0 = vec2f_1 - vec2f_2;
    vec2f_0 -= vec2f_1 - vec2f_2;

    vec2f_0 = vec2f_1 * vec2f_2;
    vec2f_0 *= vec2f_1 * vec2f_2;

    vec2f_0 = vec2f_1 / vec2f_2;
    vec2f_0 /= vec2f_1 / vec2f_2;

    dot = vec2f_0.dot(vec2f_1);
    length2 = vec2f_0.length2();
    length = vec2f_0.length();

    Vec2<float> vec2f_3 = vec2f_0.normalized();
    vec2f_0.normalize();
    //-------------------------------------------------------------------   Vec2



    //-------------------------------------------------------------------   Mat3
    Mat3<float> m1(9, 1, 1, 1, 9, 1, 1, 1, 9);
    m1.invert();
    m1.inverse();

    Mat4<float> m2(9, 1, 1, 1, 1, 9, 1, 1, 1, 1, 9, 1, 1, 1, 1, 9);
    m2.invert();
    m2.inverse();
    //-------------------------------------------------------------------   Mat3

}