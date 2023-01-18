#ifndef RED_MATH_H
#define RED_MATH_H

#include <cmath>
#include <concepts>

template<typename T>
concept Arithmetic = std::integral<T> || std::floating_point<T>;

#pragma region Forward declaration
template<Arithmetic T>
struct Vec4;

template<Arithmetic T>
struct Vec3;

template<Arithmetic T>
struct Vec2;

template<Arithmetic T, bool IsColumnMajor = true>
struct Mat3;

template<Arithmetic T, bool IsColumnMajor = true>
struct Mat4;
#pragma endregion

#pragma region Vec2
template<Arithmetic T>
struct Vec2
{
	union { T x, r, s; };
	union { T y, g, t; };

#pragma region Constructors
	Vec2() = default;
	explicit Vec2(T value)
		: x(value), y(value)
	{}
	Vec2(T x, T y)
		: x(x), y(y)
	{}
#pragma endregion

#pragma region Operators
	Vec2 operator+(const Vec2& other) const
	{
		return {
			x + other.x,
			y + other.y
		};
	}

	Vec2& operator+=(const Vec2& other)
	{
		x += other.x;
		y += other.y;

		return *this;
	}

	template <Arithmetic K>
	Vec2 operator+(K scalar) const
	{
		return {
			x + scalar,
			y + scalar
		};
	}

	template <Arithmetic K>
	Vec2& operator+=(K scalar)
	{
		x += scalar;
		y += scalar;

		return *this;
	}

	Vec2 operator-(const Vec2& other) const
	{
		return {
			x - other.x,
			y - other.y
		};
	}

	Vec2& operator-=(const Vec2& other)
	{
		x -= other.x;
		y -= other.y;

		return *this;
	}

	template <Arithmetic K>
	Vec2 operator-(K scalar) const
	{
		return {
			x - scalar,
			y - scalar
		};
	}

	template <Arithmetic K>
	Vec2& operator-=(K scalar)
	{
		x -= scalar;
		y -= scalar;

		return *this;
	}

	Vec2 operator*(const Vec2& other) const
	{
		return {
			x * other.x,
			y * other.y
		};
	}

	Vec2& operator*=(const Vec2& other)
	{
		x *= other.x;
		y *= other.y;

		return *this;
	}

	template <Arithmetic K>
	Vec2 operator*(K scalar) const
	{
		return {
			x * scalar,
			y * scalar
		};
	}

	template <Arithmetic K>
	Vec2& operator*=(K scalar)
	{
		x *= scalar;
		y *= scalar;

		return *this;
	}

	Vec2 operator/(const Vec2& other) const
	{
		return {
			x / other.x,
			y / other.y
		};
	}

	Vec2& operator/=(const Vec2& other)
	{
		x /= other.x;
		y /= other.y;

		return *this;
	}

	template <Arithmetic K>
	Vec2 operator/(K scalar) const
	{
		return {
			x / scalar,
			y / scalar
		};
	}

	template <Arithmetic K>
	Vec2& operator/=(K scalar)
	{
		x /= scalar;
		y /= scalar;

		return *this;
	}

	T& operator[](int index)
	{
		return *(&x + index);
	}

	const T& operator[](int index) const
	{
		return *(&x + index);
	}
#pragma endregion

	T dot(const Vec2& other) const
	{
		return x * other.x + y * other.y;
	}

	T length2() const
	{
		return x * x + y * y;
	}

	template<Arithmetic K = float>
	K length() const
	{
		if constexpr (std::is_same_v<K, float>)
		{
			return std::sqrtf(length2());
		}

		if constexpr (std::is_same_v<K, long double>)
		{
			return std::sqrtl(length2());
		}

		return std::sqrt(length2());
	}

	T distance2(const Vec2& other)
	{
		return (*this - other).length2();
	}

	template<Arithmetic K = float>
	K distance(const Vec2& other)
	{
		return (*this - other).length<K>();
	}

	void normalize()
	{
		if constexpr (std::is_same_v<T, double>)
		{
			double _length = length<double>();

			x /= _length;
			y /= _length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double _length = length<long double>();

			x /= _length;
			y /= _length;
		}
		else
		{
			float _length = length<float>();

			x /= _length;
			y /= _length;
		}
	}

	[[nodiscard]] Vec2<T> normalized() const
	{
		Vec2<T> v = *this;
		if constexpr (std::is_same_v<T, double>)
		{
			double length = v.template length<double>();

			v.x /= length;
			v.y /= length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double length = v.template length<long double>();

			v.x /= length;
			v.y /= length;
		}
		else
		{
			float length = v.template length<float>();

			v.x /= length;
			v.y /= length;
		}

		return v;
	}
};
#pragma endregion

#pragma region Vec3
template<Arithmetic T>
struct Vec3
{
	union { T x, r, s; };
	union { T y, g, t; };
	union { T z, b, p; };

#pragma region Constructors
	Vec3() = default;
	explicit Vec3(T value)
		: x(value), y(value), z(value)
	{}
	Vec3(T x, T y)
		: x(x), y(y), z{}
	{}
	Vec3(T x, T y, T z)
		: x(x), y(y), z(z)
	{}
	template<Arithmetic K>
	Vec3(const Vec2<K>& v);
#pragma endregion

#pragma region Operators
	Vec3 operator+(const Vec3& other) const
	{
		return {
			x + other.x,
			y + other.y,
			z + other.z
		};
	}

	Vec3& operator+=(const Vec3& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;

		return *this;
	}

	template <Arithmetic K>
	Vec3 operator+(K scalar) const
	{
		return {
			x + scalar,
			y + scalar,
			z + scalar
		};
	}

	template <Arithmetic K>
	Vec3& operator+=(K scalar)
	{
		x += scalar;
		y += scalar;
		z += scalar;

		return *this;
	}

	Vec3 operator-(const Vec3& other) const
	{
		return {
			x - other.x,
			y - other.y,
			z - other.z
		};
	}

	Vec3& operator-=(const Vec3& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;

		return *this;
	}

	template <Arithmetic K>
	Vec3 operator-(K scalar) const
	{
		return {
			x - scalar,
			y - scalar,
			z - scalar
		};
	}

	template <Arithmetic K>
	Vec3& operator-=(K scalar)
	{
		x -= scalar;
		y -= scalar;
		z -= scalar;

		return *this;
	}

	Vec3 operator*(const Vec3& other) const
	{
		return {
			x * other.x,
			y * other.y,
			z * other.z
		};
	}

	Vec3& operator*=(const Vec3& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;

		return *this;
	}

	template <Arithmetic K>
	Vec3 operator*(K scalar) const
	{
		return {
			x * scalar,
			y * scalar,
			z * scalar
		};
	}

	template <Arithmetic K>
	Vec3& operator*=(K scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;

		return *this;
	}

	Vec3 operator/(const Vec3& other) const
	{
		return {
			x / other.x,
			y / other.y,
			z / other.z
		};
	}

	Vec3& operator/=(const Vec3& other)
	{
		x /= other.x;
		y /= other.y;
		z /= other.z;

		return *this;
	}

	template <Arithmetic K>
	Vec3 operator/(K scalar) const
	{
		return {
			x / scalar,
			y / scalar,
			z / scalar
		};
	}

	template <Arithmetic K>
	Vec3& operator/=(K scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;

		return *this;
	}

	T& operator[](int index)
	{
		return *(&x + index);
	}

	const T& operator[](int index) const
	{
		return *(&x + index);
	}
#pragma endregion

	T dot(const Vec3& other) const
	{
		return x * other.x + y * other.y + z * other.z;
	}

	T length2() const
	{
		return x * x + y * y + z * z;
	}

	template<Arithmetic K = float>
	K length() const
	{
		if constexpr (std::is_same_v<K, float>)
		{
			return std::sqrtf(length2());
		}

		if constexpr (std::is_same_v<K, long double>)
		{
			return std::sqrtl(length2());
		}

		return std::sqrt(length2());
	}

	T distance2(const Vec3& other)
	{
		return (*this - other).length2();
	}

	template<Arithmetic K = float>
	K distance(const Vec3& other)
	{
		return (*this - other).length<K>();
	}

	void normalize()
	{
		if constexpr (std::is_same_v<T, double>)
		{
			double _length = length<double>();

			x /= _length;
			y /= _length;
			z /= _length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double _length = length<long double>();

			x /= _length;
			y /= _length;
			z /= _length;
		}
		else
		{
			float _length = length<float>();

			x /= _length;
			y /= _length;
			z /= _length;
		}
	}

	[[nodiscard]] Vec3 normalized() const
	{
		Vec3 v = *this;
		if constexpr (std::is_same_v<T, double>)
		{
			double length = v.template length<double>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double length = v.template length<long double>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
		}
		else
		{
			float length = v.template length<float>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
		}

		return v;
	}
	template<typename K>
	Vec3 cross(const Vec3<K>& other) const
	{
		return {
			y * other.z - z * other.y,
			z * other.x - x * other.z,
			x * other.y - y * other.x
		};
	}
};
#pragma endregion

#pragma region Vec4
template<Arithmetic T>
struct Vec4
{
	union { T x, r, s; };
	union { T y, g, t; };
	union { T z, b, p; };
	union { T w, a, q; };

#pragma region Constructors
	Vec4() = default;
	explicit Vec4(T value)
		: x(value), y(value), z(value), w(value)
	{}
	Vec4(T x, T y)
		: x(x), y(y), z{}, w{}
	{}
	Vec4(T x, T y, T z)
		: x(x), y(y), z(z), w{}
	{}
	Vec4(T x, T y, T z, T w)
		: x(x), y(y), z(z), w(w)
	{}

	template<Arithmetic K>
	Vec4(const Vec2<K>& v);

	template<Arithmetic K>
	Vec4(const Vec3<K>& v);
#pragma endregion

#pragma region Operators
	Vec4 operator+(const Vec4& other) const
	{
		return {
			x + other.x,
			y + other.y,
			z + other.z,
			w + other.w
		};
	}

	Vec4& operator+=(const Vec4& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;

		return *this;
	}

	template <Arithmetic K>
	Vec4 operator+(K scalar) const
	{
		return {
			x + scalar,
			y + scalar,
			z + scalar,
			w + scalar
		};
	}

	template <Arithmetic K>
	Vec4& operator+=(K scalar)
	{
		x += scalar;
		y += scalar;
		z += scalar;
		w += scalar;

		return *this;
	}

	Vec4 operator-(const Vec4& other) const
	{
		return {
			x - other.x,
			y - other.y,
			z - other.z,
			w - other.w
		};
	}

	Vec4& operator-=(const Vec4& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;

		return *this;
	}

	template <Arithmetic K>
	Vec4 operator-(K scalar) const
	{
		return {
			x - scalar,
			y - scalar,
			z - scalar,
			w - scalar
		};
	}

	template <Arithmetic K>
	Vec4& operator-=(K scalar)
	{
		x -= scalar;
		y -= scalar;
		z -= scalar;
		w -= scalar;

		return *this;
	}

	Vec4 operator*(const Vec4& other) const
	{
		return {
			x * other.x,
			y * other.y,
			z * other.z,
			w * other.w
		};
	}

	Vec4& operator*=(const Vec4& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;

		return *this;
	}

	template <Arithmetic K>
	Vec4 operator*(K scalar) const
	{
		return {
			x * scalar,
			y * scalar,
			z * scalar,
			w * scalar
		};
	}

	template <Arithmetic K>
	Vec4& operator*=(K scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;

		return *this;
	}

	Vec4 operator/(const Vec4& other) const
	{
		return {
			x / other.x,
			y / other.y,
			z / other.z,
			w / other.w
		};
	}

	Vec4& operator/=(const Vec4& other)
	{
		x /= other.x;
		y /= other.y;
		z /= other.z;
		w /= other.w;

		return *this;
	}

	template <Arithmetic K>
	Vec4 operator/(K scalar) const
	{
		return {
			x / scalar,
			y / scalar,
			z / scalar,
			w / scalar
		};
	}

	template <Arithmetic K>
	Vec4& operator/=(K scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;

		return *this;
	}

	T& operator[](int index)
	{
		return *(&x + index);
	}

	const T& operator[](int index) const
	{
		return *(&x + index);
	}
#pragma endregion

	T dot(const Vec4& other) const
	{
		return x * other.x + y * other.y + z * other.z + w * other.w;
	}

	T length2() const
	{
		return x * x + y * y + z * z + w * w;
	}

	template<Arithmetic K = float>
	K length() const
	{
		if constexpr (std::is_same_v<K, float>)
		{
			return std::sqrtf(length2());
		}

		if constexpr (std::is_same_v<K, long double>)
		{
			return std::sqrtl(length2());
		}

		return std::sqrt(length2());
	}

	T distance2(const Vec4& other)
	{
		return (*this - other).length2();
	}

	template<Arithmetic K = float>
	K distance(const Vec4& other)
	{
		return (*this - other).template length<K>();
	}

	void normalize()
	{
		if constexpr (std::is_same_v<T, double>)
		{
			double _length = length<double>();

			x /= _length;
			y /= _length;
			z /= _length;
			w /= _length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double _length = length<long double>();

			x /= _length;
			y /= _length;
			z /= _length;
			w /= _length;
		}
		else
		{
			float _length = length<float>();

			x /= _length;
			y /= _length;
			z /= _length;
			w /= _length;
		}
	}

	[[nodiscard]] Vec4<T> normalized() const
	{
		Vec4<T> v = *this;
		if constexpr (std::is_same_v<T, double>)
		{
			double length = v.template length<double>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
			v.w /= length;
		}
		else if constexpr (std::is_same_v<T, long double>)
		{
			long double length = v.template length<long double>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
			v.w /= length;
		}
		else
		{
			float length = v.template length<float>();

			v.x /= length;
			v.y /= length;
			v.z /= length;
			v.w /= length;
		}

		return v;
	}
};
#pragma endregion

#pragma region Mat3
template<Arithmetic T, bool IsColumnMajor>
struct Mat3
{
	Vec3<T> m0, m1, m2;

#pragma region Constructors
	Mat3() = default;

	explicit Mat3(T value)
		: m0(value), m1(value), m2(value)
	{}

	Mat3(const Vec3<T>& m0, const Vec3<T>& m1, const Vec3<T>& m2)
		: m0(m0), m1(m1), m2(m2)
	{}

	Mat3(T m00, T m01, T m02,
		T m10, T m11, T m12,
		T m20, T m21, T m22)
	{
		if constexpr (IsColumnMajor)
		{
			m0 = { m00, m10, m20 };
			m1 = { m01, m11, m21 };
			m2 = { m02, m12, m22 };
		}
		else
		{
			m0 = { m00, m01, m02 };
			m1 = { m10, m11, m12 };
			m2 = { m20, m21, m22 };
		}
	}

	static Mat3 identity()
	{
		return {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1
		};
	}
#pragma endregion

#pragma region Operators
	Mat3 operator+(const Mat3& other) const
	{
		return {
			m0 + other.m0,
			m1 + other.m1,
			m2 + other.m2
		};
	}

	Mat3& operator+=(const Mat3& other)
	{
		m0 += other.m0;
		m1 += other.m1;
		m2 += other.m2;

		return *this;
	}

	Mat3 operator-(const Mat3& other) const
	{
		return {
			m0 - other.m0,
			m1 - other.m1,
			m2 - other.m2
		};
	}

	Mat3& operator-=(const Mat3& other)
	{
		m0 -= other.m0;
		m1 -= other.m1;
		m2 -= other.m2;

		return *this;
	}

	Mat3 operator*(const Mat3& other) const
	{
		if constexpr (IsColumnMajor)
		{
			return
			{
				{
					(*this)(0, 0) * other(0, 0) + (*this)(0, 1) * other(1, 0) + (*this)(0, 2) * other(2, 0),
					(*this)(1, 0) * other(0, 0) + (*this)(1, 1) * other(1, 0) + (*this)(1, 2) * other(2, 0),
					(*this)(2, 0) * other(0, 0) + (*this)(2, 1) * other(1, 0) + (*this)(2, 2) * other(2, 0)
				},

				{
					(*this)(0, 0)* other(0, 1) + (*this)(0, 1) * other(1, 1) + (*this)(0, 2) * other(2, 1),
					(*this)(1, 0)* other(0, 1) + (*this)(1, 1) * other(1, 1) + (*this)(1, 2) * other(2, 1),
					(*this)(2, 0)* other(0, 1) + (*this)(2, 1) * other(1, 1) + (*this)(2, 2) * other(2, 1)
				},

				{
					(*this)(0, 0) * other(0, 2) + (*this)(0, 1) * other(1, 2) + (*this)(0, 2) * other(2, 2),
					(*this)(1, 0) * other(0, 2) + (*this)(1, 1) * other(1, 2) + (*this)(1, 2) * other(2, 2),
					(*this)(2, 0) * other(0, 2) + (*this)(2, 1) * other(1, 2) + (*this)(2, 2) * other(2, 2)
				}
			};
		}
		else
		{
			return
			{
				{
					(*this)(0, 0) * other(0, 0) + (*this)(0, 1) * other(1, 0) + (*this)(0, 2) * other(2, 0),
					(*this)(0, 0)* other(0, 1) + (*this)(0, 1) * other(1, 1) + (*this)(0, 2) * other(2, 1),
					(*this)(0, 0)* other(0, 2) + (*this)(0, 1) * other(1, 2) + (*this)(0, 2) * other(2, 2)
				},

				{
					(*this)(1, 0) * other(0, 0) + (*this)(1, 1) * other(1, 0) + (*this)(1, 2) * other(2, 0),
					(*this)(1, 0)* other(0, 1) + (*this)(1, 1) * other(1, 1) + (*this)(1, 2) * other(2, 1),
					(*this)(1, 0)* other(0, 2) + (*this)(1, 1) * other(1, 2) + (*this)(1, 2) * other(2, 2)
				},

				{
					(*this)(2, 0) * other(0, 0) + (*this)(2, 1) * other(1, 0) + (*this)(2, 2) * other(2, 0),
					(*this)(2, 0)* other(0, 1) + (*this)(2, 1) * other(1, 1) + (*this)(2, 2) * other(2, 1),
					(*this)(2, 0)* other(0, 2) + (*this)(2, 1) * other(1, 2) + (*this)(2, 2) * other(2, 2)
				}
			};
		}
	}

	Mat3& operator*=(const Mat3& other)
	{
		return *this = *this * other;
		//if constexpr (IsColumnMajor)
		//{
		//	m0 =
		//	{
		//		(*this)(0, 0] * other(0, 0] + (*this)(0, 1] * other(1, 0] + (*this)(0, 2] * other(2, 0],
		//		(*this)(1, 0] * other(0, 0] + (*this)(1, 1] * other(1, 0] + (*this)(1, 2] * other(2, 0],
		//		(*this)(2, 0] * other(0, 0] + (*this)(2, 1] * other(1, 0] + (*this)(2, 2] * other(2, 0],
		//	};
		//
		//	m1 =
		//	{
		//		(*this)(0, 0] * other(0, 1] + (*this)(0, 1] * other(1, 1] + (*this)(0, 2] * other(2, 1],
		//		(*this)(1, 0] * other(0, 1] + (*this)(1, 1] * other(1, 1] + (*this)(1, 2] * other(2, 1],
		//		(*this)(2, 0] * other(0, 1] + (*this)(2, 1] * other(1, 1] + (*this)(2, 2] * other(2, 1],
		//	};
		//
		//	m2 =
		//	{
		//		(*this)(0, 0] * other(0, 2] + (*this)(0, 1] * other(1, 2] + (*this)(0, 2] * other(2, 2],
		//		(*this)(1, 0] * other(0, 2] + (*this)(1, 1] * other(1, 2] + (*this)(1, 2] * other(2, 2],
		//		(*this)(2, 0] * other(0, 2] + (*this)(2, 1] * other(1, 2] + (*this)(2, 2] * other(2, 2],
		//	};
		//
		//}
		//else
		//{
		//
		//	m0 = {
		//			(*this)[0, 0] * other[0, 0] + (*this)[0, 1] * other[1, 0] + (*this)[0, 2] * other[2, 0],
		//			(*this)[0, 0] * other[0, 1] + (*this)[0, 1] * other[1, 1] + (*this)[0, 2] * other[2, 1],
		//			(*this)[0, 0] * other[0, 2] + (*this)[0, 1] * other[1, 2] + (*this)[0, 2] * other[2, 2],
		//	};
		//
		//	m1 = {
		//			(*this)[1, 0] * other[0, 0] + (*this)[1, 1] * other[1, 0] + (*this)[1, 2] * other[2, 0],
		//			(*this)[1, 0] * other[0, 1] + (*this)[1, 1] * other[1, 1] + (*this)[1, 2] * other[2, 1],
		//			(*this)[1, 0] * other[0, 2] + (*this)[1, 1] * other[1, 2] + (*this)[1, 2] * other[2, 2],
		//	};
		//
		//	m2 = {
		//			(*this)[2, 0] * other[0, 0] + (*this)[2, 1] * other[1, 0] + (*this)[2, 2] * other[2, 0],
		//			(*this)[2, 0] * other[0, 1] + (*this)[2, 1] * other[1, 1] + (*this)[2, 2] * other[2, 1],
		//			(*this)[2, 0] * other[0, 2] + (*this)[2, 1] * other[1, 2] + (*this)[2, 2] * other[2, 2]
		//	};
		//}
		//
		//return *this;
	}

	template<Arithmetic K>
	Mat3 operator*(K scalar) const
	{
		return *this * Mat3(scalar, 0, 0,
							0, scalar, 0,
							0, 0, scalar);
	}

	template<Arithmetic K>
	Mat3 operator*=(K scalar)
	{
		*this = *this * Mat3(scalar, 0, 0,
							0, scalar, 0,
							0, 0, scalar);
		return *this;
	}

	T& operator()(int row, int column)
	{
		if constexpr (IsColumnMajor)
		{
			return *(&m0[0] + 3 * column + row);
		}
		else
		{
			return *(&m0[0] + 3 * row + column);
		}
	}

	const T& operator()(int row, int column) const
	{
		if constexpr (IsColumnMajor)
		{
			return *(&m0[0] + 3 * column + row);
		}
		else
		{
			return *(&m0[0] + 3 * row + column);
		}
	}
#pragma endregion

	Mat3 transposed() const
	{
		return {
			(*this)(0, 0), (*this)(1, 0), (*this)(2, 0),
			(*this)(0, 1), (*this)(1, 1), (*this)(2, 1),
			(*this)(0, 2), (*this)(1, 2), (*this)(2, 2)
		};
	}

	void transpose()
	{
		*this = transposed();
	}

	T determinant() const
	{
		T m00 = (*this)(0, 0);
		T m01 = (*this)(0, 1);
		T m02 = (*this)(0, 2);

		T m10 = (*this)(1, 0);
		T m11 = (*this)(1, 1);
		T m12 = (*this)(1, 2);

		T m20 = (*this)(2, 0);
		T m21 = (*this)(2, 1);
		T m22 = (*this)(2, 2);

		return m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
	}

	Mat3 inverse() const
	{
		T m00 = (*this)(0, 0);
		T m01 = (*this)(0, 1);
		T m02 = (*this)(0, 2);
					   
		T m10 = (*this)(1, 0);
		T m11 = (*this)(1, 1);
		T m12 = (*this)(1, 2);
					   
		T m20 = (*this)(2, 0);
		T m21 = (*this)(2, 1);
		T m22 = (*this)(2, 2);
		
		T min00 = m11 * m22 - m12 * m21;
		T min01 = m10 * m22 - m12 * m20;
		T min02 = m10 * m21 - m11 * m20;

		T min10 = m01 * m22 - m02 * m21;
		T min11 = m00 * m22 - m02 * m20;
		T min12 = m00 * m21 - m01 * m20;

		T min20 = m01 * m12 - m02 * m11;
		T min21 = m00 * m12 - m02 * m10;
		T min22 = m00 * m11 - m01 * m10;

		Mat3 minMat = {
			 min00, -min01,  min02,
			-min10,  min11, -min12,
			 min20, -min21,  min22
		};

		minMat.transpose();
		return minMat * (1.0f / determinant());
	}
	void invert()
	{
		*this = inverse();
	}

};
#pragma endregion

#pragma region Mat4
template<Arithmetic T, bool IsColumnMajor>
struct Mat4
{
	Vec4<T> m0, m1, m2, m3;

#pragma region Constructors
	Mat4() = default;

	explicit Mat4(T value)
		: m0(value), m1(value), m2(value), m3(value)
	{}

	Mat4(const Vec4<T>& m0, const Vec4<T>& m1, const Vec4<T>& m2, const Vec4<T>& m3)
		: m0(m0), m1(m1), m2(m2), m3(m3)
	{}

	Mat4(T m00, T m01, T m02, T m03,
		T m10, T m11, T m12, T m13,
		T m20, T m21, T m22, T m23,
		T m30, T m31, T m32, T m33)
	{
		if constexpr (IsColumnMajor)
		{
			m0 = { m00, m10, m20, m30 };
			m1 = { m01, m11, m21, m31 };
			m2 = { m02, m12, m22, m32 };
			m3 = { m03, m13, m23, m33 };
		}
		else
		{
			m0 = { m00, m01, m02, m03 };
			m1 = { m10, m11, m12, m13 };
			m2 = { m20, m21, m22, m23 };
			m3 = { m30, m31, m32, m33 };
		}
	}

#pragma endregion

#pragma region Operators
	Mat4 operator+(const Mat4& other) const
	{
		return {
			m0 + other.m0,
			m1 + other.m1,
			m2 + other.m2,
			m3 + other.m3
		};
	}

	Mat4& operator+=(const Mat4& other)
	{
		m0 += other.m0;
		m1 += other.m1;
		m2 += other.m2;
		m3 += other.m3;

		return *this;
	}

	Mat4 operator-(const Mat4& other) const
	{
		return {
			m0 - other.m0,
			m1 - other.m1,
			m2 - other.m2,
			m3 - other.m3
		};
	}

	Mat4& operator-=(const Mat4& other)
	{
		m0 -= other.m0;
		m1 -= other.m1;
		m2 -= other.m2;
		m3 -= other.m3;

		return *this;
	}

	Mat4 operator*(const Mat4& other) const
	{
		if constexpr (IsColumnMajor)
		{
			return
			{
				{
					(*this)(0, 0) * other(0, 0) + (*this)(0, 1) * other(1, 0) + (*this)(0, 2) * other(2, 0) + (*this)(0, 3) * other(3, 0),
					(*this)(1, 0) * other(0, 0) + (*this)(1, 1) * other(1, 0) + (*this)(1, 2) * other(2, 0) + (*this)(1, 3) * other(3, 0),
					(*this)(2, 0) * other(0, 0) + (*this)(2, 1) * other(1, 0) + (*this)(2, 2) * other(2, 0) + (*this)(2, 3) * other(3, 0),
					(*this)(3, 0) * other(0, 0) + (*this)(3, 1) * other(1, 0) + (*this)(3, 2) * other(2, 0) + (*this)(3, 3) * other(3, 0)
				},	
					
				{	
					(*this)(0, 0) * other(0, 1) + (*this)(0, 1) * other(1, 1) + (*this)(0, 2) * other(2, 1) + (*this)(0, 3) * other(3, 1),
					(*this)(1, 0) * other(0, 1) + (*this)(1, 1) * other(1, 1) + (*this)(1, 2) * other(2, 1) + (*this)(1, 3) * other(3, 1),
					(*this)(2, 0) * other(0, 1) + (*this)(2, 1) * other(1, 1) + (*this)(2, 2) * other(2, 1) + (*this)(2, 3) * other(3, 1),
					(*this)(3, 0) * other(0, 1) + (*this)(3, 1) * other(1, 1) + (*this)(3, 2) * other(2, 1) + (*this)(3, 3) * other(3, 1)
				},	
					
				{	
					(*this)(0, 0) * other(0, 2) + (*this)(0, 1) * other(1, 2) + (*this)(0, 2) * other(2, 2) + (*this)(0, 3) * other(3, 2),
					(*this)(1, 0) * other(0, 2) + (*this)(1, 1) * other(1, 2) + (*this)(1, 2) * other(2, 2) + (*this)(1, 3) * other(3, 2),
					(*this)(2, 0) * other(0, 2) + (*this)(2, 1) * other(1, 2) + (*this)(2, 2) * other(2, 2) + (*this)(2, 3) * other(3, 2),
					(*this)(3, 0) * other(0, 2) + (*this)(3, 1) * other(1, 2) + (*this)(3, 2) * other(2, 2) + (*this)(3, 3) * other(3, 2)
				},		 				 
						 				 
				{		 				 
					(*this)(0, 0) * other(0, 3) + (*this)(0, 1) * other(1, 3) + (*this)(0, 2) * other(2, 3) + (*this)(0, 3) * other(3, 3),
					(*this)(1, 0) * other(0, 3) + (*this)(1, 1) * other(1, 3) + (*this)(1, 2) * other(2, 3) + (*this)(1, 3) * other(3, 3),
					(*this)(2, 0) * other(0, 3) + (*this)(2, 1) * other(1, 3) + (*this)(2, 2) * other(2, 3) + (*this)(2, 3) * other(3, 3),
					(*this)(3, 0) * other(0, 3) + (*this)(3, 1) * other(1, 3) + (*this)(3, 2) * other(2, 3) + (*this)(3, 3) * other(3, 3)
				}						 
			};							 
		}								 
		else							 
		{								 
			return						 
			{							 
				{						 
					(*this)(0, 0) * other(0, 0) + (*this)(0, 1) * other(1, 0) + (*this)(0, 2) * other(2, 0) + (*this)(0, 3) * other(3, 0),
					(*this)(0, 0) * other(0, 1) + (*this)(0, 1) * other(1, 1) + (*this)(0, 2) * other(2, 1) + (*this)(0, 3) * other(3, 1),
					(*this)(0, 0) * other(0, 2) + (*this)(0, 1) * other(1, 2) + (*this)(0, 2) * other(2, 2) + (*this)(0, 3) * other(3, 2),
					(*this)(0, 0) * other(0, 3) + (*this)(0, 1) * other(1, 3) + (*this)(0, 2) * other(2, 3) + (*this)(0, 3) * other(3, 3)
				},		  				  
						  				  
				{		  				  
					(*this)(1, 0) * other(0, 0) + (*this)(1, 1) * other(1, 0) + (*this)(1, 2) * other(2, 0) + (*this)(1, 3) * other(3, 0),
					(*this)(1, 0) * other(0, 1) + (*this)(1, 1) * other(1, 1) + (*this)(1, 2) * other(2, 1) + (*this)(1, 3) * other(3, 1),
					(*this)(1, 0) * other(0, 2) + (*this)(1, 1) * other(1, 2) + (*this)(1, 2) * other(2, 2) + (*this)(1, 3) * other(3, 2),
					(*this)(1, 0) * other(0, 3) + (*this)(1, 1) * other(1, 3) + (*this)(1, 2) * other(2, 3) + (*this)(1, 3) * other(3, 3)
				},						  
										  
				{						  
					(*this)(2, 0) * other(0, 0) + (*this)(2, 1) * other(1, 0) + (*this)(2, 2) * other(2, 0) + (*this)(2, 3) * other(3, 0),
					(*this)(2, 0) * other(0, 1) + (*this)(2, 1) * other(1, 1) + (*this)(2, 2) * other(2, 1) + (*this)(2, 3) * other(3, 1),
					(*this)(2, 0) * other(0, 2) + (*this)(2, 1) * other(1, 2) + (*this)(2, 2) * other(2, 2) + (*this)(2, 3) * other(3, 2),
					(*this)(2, 0) * other(0, 3) + (*this)(2, 1) * other(1, 3) + (*this)(2, 2) * other(2, 3) + (*this)(2, 3) * other(3, 3)
				},		 				   
						 				   
				{		 				   
					(*this)(3, 0) * other(0, 0) + (*this)(3, 1) * other(1, 0) + (*this)(3, 2) * other(2, 0) + (*this)(3, 3) * other(3, 0),
					(*this)(3, 0) * other(0, 1) + (*this)(3, 1) * other(1, 1) + (*this)(3, 2) * other(2, 1) + (*this)(3, 3) * other(3, 1),
					(*this)(3, 0) * other(0, 2) + (*this)(3, 1) * other(1, 2) + (*this)(3, 2) * other(2, 2) + (*this)(3, 3) * other(3, 2),
					(*this)(3, 0) * other(0, 3) + (*this)(3, 1) * other(1, 3) + (*this)(3, 2) * other(2, 3) + (*this)(3, 3) * other(3, 3)
				}
			};
		}
	}

	Mat4& operator*=(const Mat4& other)
	{
		return *this = *this * other;
	}

	template<Arithmetic K>
	Mat4 operator*(K scalar) const
	{
		return *this * Mat4(scalar, 0, 0, 0,
							0, scalar, 0, 0,
							0, 0, scalar, 0,
							0, 0, 0, scalar);
	}

	template<Arithmetic K>
	Mat4 operator*=(K scalar)
	{
		*this = *this * Mat4(scalar, 0, 0, 0,
							0, scalar, 0, 0,
							0, 0, scalar, 0,
							0, 0, 0, scalar);
		return *this;
	}

	T& operator()(int row, int column)
	{
		if constexpr (IsColumnMajor)
		{
			return *(&m0[0] + 4 * column + row);
		}
		else
		{
			return *(&m0[0] + 4 * row + column);
		}
	}

	const T& operator()(int row, int column) const
	{
		if constexpr (IsColumnMajor)
		{
			return *(&m0[0] + 4 * column + row);
		}
		else
		{
			return *(&m0[0] + 4 * row + column);
		}
	}
#pragma endregion

	Mat4 transposed() const
	{
		return {
			(*this)(0, 0), (*this)(1, 0), (*this)(2, 0), (*this)(3, 0),
			(*this)(0, 1), (*this)(1, 1), (*this)(2, 1), (*this)(3, 1),
			(*this)(0, 2), (*this)(1, 2), (*this)(2, 2), (*this)(3, 2),
			(*this)(0, 3), (*this)(1, 3), (*this)(2, 3), (*this)(3, 3)
		};
	}

	void transpose()
	{
		*this = transposed();
	}

	T determinant() const
	{
		T a = (*this)(0, 0);
		T b = (*this)(0, 1);
		T c = (*this)(0, 2);
		T d = (*this)(0, 3);
		T e = (*this)(1, 0);
		T f = (*this)(1, 1);
		T g = (*this)(1, 2);
		T h = (*this)(1, 3);
		T i = (*this)(2, 0);
		T j = (*this)(2, 1);
		T k = (*this)(2, 2);
		T l = (*this)(2, 3);
		T m = (*this)(3, 0);
		T n = (*this)(3, 1);
		T o = (*this)(3, 2);
		T p = (*this)(3, 3);

		auto aa = Mat3(f, j, n, g, k, o, h, l, p).determinant();
		auto bb = Mat3(e, i, m, g, k, o, h, l, p).determinant();
		auto cc = Mat3(e, i, m, f, j, n, h, l, p).determinant();
		auto dd = Mat3(e, i, m, f, j, n, g, k, o).determinant();
		return aa * a - bb * b + cc * c - dd * d;
	}
	[[nodiscard]] Mat4 inverse() const
	{
		T m00 = (*this)(0, 0);
		T m01 = (*this)(0, 1);
		T m02 = (*this)(0, 2);
		T m03 = (*this)(0, 3);

		T m10 = (*this)(1, 0);
		T m11 = (*this)(1, 1);
		T m12 = (*this)(1, 2);
		T m13 = (*this)(1, 3);

		T m20 = (*this)(2, 0);
		T m21 = (*this)(2, 1);
		T m22 = (*this)(2, 2);
		T m23 = (*this)(2, 3);

		T m30 = (*this)(3, 0);
		T m31 = (*this)(3, 1);
		T m32 = (*this)(3, 2);
		T m33 = (*this)(3, 3);

		T min00 = Mat3<T, false>(m11, m12, m13, m21, m22, m23, m31, m32, m33).determinant();
		T min01 = Mat3<T, false>(m10, m12, m13, m20, m22, m23, m30, m32, m33).determinant();
		T min02 = Mat3<T, false>(m10, m11, m13, m20, m21, m23, m30, m31, m33).determinant();
		T min03 = Mat3<T, false>(m10, m11, m12, m20, m21, m22, m30, m31, m32).determinant();

		T min10 = Mat3<T, false>(m01, m02, m03, m21, m22, m23, m31, m32, m33).determinant();
		T min11 = Mat3<T, false>(m00, m02, m03, m20, m22, m23, m30, m32, m33).determinant();
		T min12 = Mat3<T, false>(m00, m01, m03, m20, m21, m23, m30, m31, m33).determinant();
		T min13 = Mat3<T, false>(m00, m01, m02, m20, m21, m22, m30, m31, m32).determinant();

		T min20 = Mat3<T, false>(m01, m02, m03, m11, m12, m13, m31, m32, m33).determinant();
		T min21 = Mat3<T, false>(m00, m02, m03, m10, m12, m13, m30, m32, m33).determinant();
		T min22 = Mat3<T, false>(m00, m01, m03, m10, m11, m13, m30, m31, m33).determinant();
		T min23 = Mat3<T, false>(m00, m01, m02, m10, m11, m12, m30, m31, m32).determinant();

		T min30 = Mat3<T, false>(m01, m02, m03, m11, m12, m13, m21, m22, m23).determinant();
		T min31 = Mat3<T, false>(m00, m02, m03, m10, m12, m13, m20, m22, m23).determinant();
		T min32 = Mat3<T, false>(m00, m01, m03, m10, m11, m13, m20, m21, m23).determinant();
		T min33 = Mat3<T, false>(m00, m01, m02, m10, m11, m12, m20, m21, m22).determinant();

		Mat4 minMat = {
			 min00, -min01,  min02, -min03,
			-min10,  min11, -min12,  min13,
			 min20, -min21,  min22, -min23,
			-min30,  min31, -min32,  min33
		};

		minMat.transpose();
		return minMat * (1.0f / determinant());
	}
	void invert()
	{
		*this = inverse();
	}

};
#pragma endregion

#pragma region Remaining definitions

template <Arithmetic T>
template <Arithmetic K>
Vec3<T>::Vec3(const Vec2<K>& v)
	: x(v.x), y(v.y), z(0.0f)
{}

template <Arithmetic T>
template <Arithmetic K>
Vec4<T>::Vec4(const Vec2<K>& v)
	: x(v.x), y(v.y), z{}, w{}
{}

template <Arithmetic T>
template <Arithmetic K>
Vec4<T>::Vec4(const Vec3<K>& v)
	: x(v.x), y(v.y), z(v.z), w{}
{}


#pragma endregion

#endif