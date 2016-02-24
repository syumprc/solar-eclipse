/*
*
* Template Numerical Toolkit (TNT)
*
*
*/

#ifndef TNT_ARRAY1D_UTILS_H
#define TNT_ARRAY1D_UTILS_H

#include <cstdlib>
#include <cassert>

namespace TNT
{


template <class T>
std::ostream& operator<<(std::ostream &s, const Array1D<T> &A)
{
    int N=A.dim1();

#ifdef TNT_DEBUG
	s << "addr: " << (void *) &A[0] << "\n";
#endif
    s << N << "\n";
    for (int j=0; j<N; j++)
    {
       s << A[j] << "\n";
    }
    s << "\n";

    return s;
}

template <class T>
std::istream& operator>>(std::istream &s, Array1D<T> &A)
{
	int N;
	s >> N;

	Array1D<T> B(N);
	for (int i=0; i<N; i++)
		s >> B[i];
	A = B;
	return s;
}



template <class T>
Array1D<T> operator+(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] + B[i];
		}
		return C;
	}
}



template <class T>
Array1D<T> operator-(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] - B[i];
		}
		return C;
	}
}


template <class T>
Array1D<T> operator*(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] * B[i];
		}
		return C;
	}
}


template <class T>
Array1D<T> operator/(const Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() != n )
		return Array1D<T>();

	else
	{
		Array1D<T> C(n);

		for (int i=0; i<n; i++)
		{
			C[i] = A[i] / B[i];
		}
		return C;
	}
}









template <class T>
Array1D<T>&  operator+=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] += B[i];
		}
	}
	return A;
}




template <class T>
Array1D<T>&  operator-=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] -= B[i];
		}
	}
	return A;
}



template <class T>
Array1D<T>&  operator*=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] *= B[i];
		}
	}
	return A;
}




template <class T>
Array1D<T>&  operator/=(Array1D<T> &A, const Array1D<T> &B)
{
	int n = A.dim1();

	if (B.dim1() == n)
	{
		for (int i=0; i<n; i++)
		{
				A[i] /= B[i];
		}
	}
	return A;
}






} // namespace TNT

#endif
